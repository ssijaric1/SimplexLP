//  SimplexLP - Project #13: Revised Simplex in natID
//  Student: Šehzada Sijarić (19964)
//
//  VizModel.h
//  State shared between the controls and the canvas of one visualization
//  tab (2D or 3D): the user model, both solver results with recorded
//  iteration paths, the feasible-region geometry and the animation cursor.
//  Pure model code - the canvas only reads from it, the controls drive it.
#pragma once

#include "../core/LPProblem.h"
#include "../core/RevisedSimplex.h"
#include "../core/InteriorPoint.h"
#include "../core/Verify.h"
#include "../core/Polytope.h"
#include "LinearExprParser.h"

namespace lp
{

class VizModel
{
public:
    int dim = 2; // 2 or 3 structural variables

    UserModel  user;
    StandardLP std_;
    bool       hasProblem = false;

    LPResult simplex;
    LPResult ipm;

    // geometry
    std::vector<Pt2> polygon;   // 2D feasible region
    Polytope3        poly3;     // 3D feasible polytope
    double           extent = 10.0; // data extent for auto-fit / safety box

    // animation cursor over simplex.path  (0 .. path.size()-1)
    int cursor = 0;

    std::string lastError;

    explicit VizModel(int dimension) : dim(dimension) {}

    // ------------------------------------------------------------------ //
    bool build(const std::string& objText, bool maximize,
               const std::string& consText)
    {
        lastError.clear();
        hasProblem = false;
        simplex = LPResult();
        ipm     = LPResult();
        cursor  = 0;

        LinearExprParser parser(dim);
        if (!parser.parseModel(objText, maximize, consText, user))
        {
            lastError = parser.error();
            return false;
        }

        std_ = toStandardForm(user);
        hasProblem = true;

        rebuildGeometry();
        return true;
    }

    void solveBoth()
    {
        if (!hasProblem)
            return;

        SimplexOptions so;
        so.recordPath = true;
        RevisedSimplex spx;
        simplex = spx.solve(std_, so);

        IPMOptions io;
        io.recordPath = true;
        InteriorPoint ip;
        ipm = ip.solve(std_, io);

        cursor = 0;
        rebuildGeometry(); // extent may grow to include solution points
    }

    // ------------------------------------------------------------------ //
    bool hasPath() const { return !simplex.path.empty(); }
    int  pathLen() const { return (int)simplex.path.size(); }

    void resetCursor()  { cursor = 0; }
    bool stepCursor()
    {
        if (cursor + 1 < pathLen())
        {
            ++cursor;
            return true;
        }
        return false;
    }
    bool cursorAtEnd() const { return cursor + 1 >= pathLen(); }

    const IterSnapshot* currentSnap() const
    {
        if (cursor >= 0 && cursor < pathLen())
            return &simplex.path[(size_t)cursor];
        return nullptr;
    }

    // user-space point of snapshot k (first dim coordinates)
    Pt2 snap2(int k) const
    {
        const auto& s = simplex.path[(size_t)k];
        return { s.x.size() > 0 ? s.x[0] : 0.0,
                 s.x.size() > 1 ? s.x[1] : 0.0 };
    }
    Pt3 snap3(int k) const
    {
        const auto& s = simplex.path[(size_t)k];
        return { s.x.size() > 0 ? s.x[0] : 0.0,
                 s.x.size() > 1 ? s.x[1] : 0.0,
                 s.x.size() > 2 ? s.x[2] : 0.0 };
    }

    Pt2 ipmSnap2(int k) const
    {
        const auto& s = ipm.path[(size_t)k];
        return { s.x.size() > 0 ? s.x[0] : 0.0,
                 s.x.size() > 1 ? s.x[1] : 0.0 };
    }
    Pt3 ipmSnap3(int k) const
    {
        const auto& s = ipm.path[(size_t)k];
        return { s.x.size() > 0 ? s.x[0] : 0.0,
                 s.x.size() > 1 ? s.x[1] : 0.0,
                 s.x.size() > 2 ? s.x[2] : 0.0 };
    }

    // status line for the log
    std::string summary() const
    {
        if (!hasProblem)
            return "no problem";
        std::string s;
        if (simplex.status != Status::NotSolved)
        {
            s += "Simplex: ";
            s += toString(simplex.status);
            if (simplex.status == Status::Optimal)
            {
                char buf[160];
                std::snprintf(buf, sizeof(buf),
                              "  obj=%.6g  iters=%d (phase I %d)  factorizations=%d  %.2f ms",
                              std_.userObjective(&simplex.x[0]),
                              simplex.stats.iterations, simplex.stats.phase1Iters,
                              simplex.stats.factorizations, simplex.stats.msTotal);
                s += buf;
            }
            else if (!simplex.message.empty())
                s += " (" + simplex.message + ")";
        }
        if (ipm.status != Status::NotSolved)
        {
            s += "\nInterior point: ";
            s += toString(ipm.status);
            if (ipm.status == Status::Optimal)
            {
                char buf[120];
                std::snprintf(buf, sizeof(buf), "  obj=%.6g  iters=%d  %.2f ms",
                              std_.userObjective(&ipm.x[0]),
                              ipm.stats.iterations, ipm.stats.msTotal);
                s += buf;
            }
            else if (!ipm.message.empty())
                s += " (" + ipm.message + ")";
        }
        if (simplex.status == Status::Optimal && ipm.status == Status::Optimal)
        {
            KKTReport k = verifyKKT(std_, simplex.x, simplex.y);
            char buf[120];
            std::snprintf(buf, sizeof(buf),
                          "\nKKT(simplex): gap=%.1e  primal=%.1e  dual=%.1e",
                          k.dualityGap, k.primalResidual, k.dualNegativity);
            s += buf;
        }
        return s;
    }

private:
    void rebuildGeometry()
    {
        // data extent: constraint intercepts and any recorded points
        double e = 1.0;
        for (const auto& r : user.rows)
        {
            for (double a : r.a)
                if (std::fabs(a) > 1e-12)
                    e = std::max(e, std::fabs(r.b / a));
        }
        for (const auto& s : simplex.path)
            for (int j = 0; j < dim && j < (int)s.x.size(); ++j)
                e = std::max(e, std::fabs(s.x[(size_t)j]));
        for (const auto& s : ipm.path)
            for (int j = 0; j < dim && j < (int)s.x.size(); ++j)
                e = std::max(e, std::fabs(s.x[(size_t)j]));
        extent = e * 1.15 + 1e-9;

        if (dim == 2)
            polygon = feasiblePolygon2D(user, extent * 4.0);
        else
            poly3 = buildPolytope3D(user, extent * 4.0);
    }
};

} // namespace lp
