//  SimplexLP - Project #13: Revised Simplex in natID
//  Student: Šehzada Sijarić (19964)
//
//  Polytope.h
//  Geometry backing the 2D/3D visualization of the feasible region:
//
//   * 2D: the region {x : a_i^T x <= b_i, x >= 0} is obtained by
//     Sutherland-Hodgman clipping of a large box against every half-plane.
//   * 3D: polytope vertices are enumerated as intersections of all triples
//     of bounding planes; each 3x3 system is solved with natID's
//     dense::Matrix::solve() (the "dense::Matrix for column operations"
//     element of the proposal). Faces are recovered per plane by sorting
//     incident vertices by angle inside the plane.
//
//  Everything here is pure geometry - no gui includes - so the whole file
//  is exercised by the standalone test suite as well.
#pragma once

#include "LPProblem.h"
#include <dense/Matrix.h>
#include <array>

namespace lp
{

// ------------------------------------------------------------------- 2D ---
struct Pt2 { double x = 0.0, y = 0.0; };

// half-plane n.x * x + n.y * y <= rhs
struct HalfPlane2
{
    double nx = 0.0, ny = 0.0, rhs = 0.0;
    inline double eval(const Pt2& p) const { return nx * p.x + ny * p.y - rhs; }
};

// Sutherland-Hodgman: clip polygon by one half-plane (keep eval <= 0).
inline std::vector<Pt2> clipPolygon(const std::vector<Pt2>& poly, const HalfPlane2& h)
{
    std::vector<Pt2> out;
    const size_t n = poly.size();
    if (n == 0)
        return out;
    out.reserve(n + 2);

    for (size_t i = 0; i < n; ++i)
    {
        const Pt2& A = poly[i];
        const Pt2& B = poly[(i + 1) % n];
        const double da = h.eval(A);
        const double db = h.eval(B);
        const bool inA = (da <= 1e-12);
        const bool inB = (db <= 1e-12);

        if (inA)
            out.push_back(A);
        if (inA != inB)
        {
            const double t = da / (da - db); // intersection parameter
            out.push_back({A.x + t * (B.x - A.x), A.y + t * (B.y - A.y)});
        }
    }
    return out;
}

// Feasible polygon of a 2-variable user model inside [0,box]^2.
// Equality rows are treated as a pair of opposite inequalities, so the
// region degenerates to a segment, which still renders correctly.
inline std::vector<Pt2> feasiblePolygon2D(const UserModel& um, double box)
{
    std::vector<Pt2> poly = { {0,0}, {box,0}, {box,box}, {0,box} };

    // x >= 0, y >= 0 are implied by the box starting at 0.
    for (const IneqConstraint& r : um.rows)
    {
        const double a0 = r.a.size() > 0 ? r.a[0] : 0.0;
        const double a1 = r.a.size() > 1 ? r.a[1] : 0.0;
        if (r.rel == Relation::LessEq || r.rel == Relation::Equal)
            poly = clipPolygon(poly, { a0,  a1,  r.b});
        if (r.rel == Relation::GreaterEq || r.rel == Relation::Equal)
            poly = clipPolygon(poly, {-a0, -a1, -r.b});
        if (poly.empty())
            break;
    }
    return poly;
}

// Intersection of two constraint lines via dense::Matrix (2x2 solve).
inline bool intersect2D(double a0, double a1, double b0,
                        double c0, double c1, double b1, Pt2& out)
{
    dense::DblMatrix M(2, 2, nullptr, true);
    auto mm = M.getManipulator();
    mm(0, 0) = a0; mm(0, 1) = a1;
    mm(1, 0) = c0; mm(1, 1) = c1;

    dense::DblMatrix B(2, 1, nullptr, true);
    auto bb = B.getFirstColumnManipulator();
    bb(0) = b0;
    bb(1) = b1;

    if (!M.solve(B))
        return false;
    auto xx = B.getManipulator();
    out.x = xx(0, 0);
    out.y = xx(1, 0);
    return std::isfinite(out.x) && std::isfinite(out.y);
}

// ------------------------------------------------------------------- 3D ---
struct Pt3
{
    double x = 0.0, y = 0.0, z = 0.0;
    Pt3() = default;
    Pt3(double X, double Y, double Z) : x(X), y(Y), z(Z) {}
    Pt3 operator-(const Pt3& o) const { return {x - o.x, y - o.y, z - o.z}; }
    Pt3 operator+(const Pt3& o) const { return {x + o.x, y + o.y, z + o.z}; }
    Pt3 operator*(double s)     const { return {x * s, y * s, z * s}; }
};

inline double dot3(const Pt3& a, const Pt3& b)  { return a.x*b.x + a.y*b.y + a.z*b.z; }
inline Pt3    cross3(const Pt3& a, const Pt3& b)
{
    return { a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x };
}
inline double norm3(const Pt3& a) { return std::sqrt(dot3(a, a)); }

// plane n^T p <= rhs
struct Plane3
{
    Pt3    n;
    double rhs = 0.0;
    int    sourceRow = -1; // index of the user constraint, -1 for axis bounds
    inline double eval(const Pt3& p) const { return dot3(n, p) - rhs; }
};

struct Face3
{
    int plane = -1;            // index into the plane list
    std::vector<int> verts;    // ordered vertex indices (CCW seen from outside)
    Pt3 normal;                // outward unit normal
    Pt3 centroid;
};

struct Polytope3
{
    std::vector<Pt3>    verts;
    std::vector<Face3>  faces;
    std::vector<Plane3> planes;
    bool empty() const { return verts.empty(); }
};

// 3x3 plane intersection through natID dense::Matrix::solve().
inline bool intersect3Planes(const Plane3& p1, const Plane3& p2, const Plane3& p3, Pt3& out)
{
    dense::DblMatrix M(3, 3, nullptr, true);
    auto mm = M.getManipulator();
    mm(0,0)=p1.n.x; mm(0,1)=p1.n.y; mm(0,2)=p1.n.z;
    mm(1,0)=p2.n.x; mm(1,1)=p2.n.y; mm(1,2)=p2.n.z;
    mm(2,0)=p3.n.x; mm(2,1)=p3.n.y; mm(2,2)=p3.n.z;

    dense::DblMatrix B(3, 1, nullptr, true);
    auto bb = B.getFirstColumnManipulator();
    bb(0) = p1.rhs; bb(1) = p2.rhs; bb(2) = p3.rhs;

    if (!M.solve(B))
        return false;
    auto xx = B.getManipulator();
    out = { xx(0,0), xx(1,0), xx(2,0) };
    return std::isfinite(out.x) && std::isfinite(out.y) && std::isfinite(out.z);
}

// Build the feasible polytope of a 3-variable user model (<=, >=, = rows,
// x,y,z >= 0, capped by a [0,box]^3 safety box so unbounded regions still
// draw something sensible).
inline Polytope3 buildPolytope3D(const UserModel& um, double box, double tol = 1e-7)
{
    Polytope3 P;

    auto addPlane = [&](double ax, double ay, double az, double rhs, int src)
    {
        Plane3 pl;
        pl.n   = {ax, ay, az};
        pl.rhs = rhs;
        pl.sourceRow = src;
        const double len = norm3(pl.n);
        if (len < 1e-14)
            return;
        pl.n   = pl.n * (1.0 / len);
        pl.rhs = pl.rhs / len;
        P.planes.push_back(pl);
    };

    for (size_t r = 0; r < um.rows.size(); ++r)
    {
        const IneqConstraint& c = um.rows[r];
        const double a0 = c.a.size() > 0 ? c.a[0] : 0.0;
        const double a1 = c.a.size() > 1 ? c.a[1] : 0.0;
        const double a2 = c.a.size() > 2 ? c.a[2] : 0.0;
        if (c.rel == Relation::LessEq || c.rel == Relation::Equal)
            addPlane( a0,  a1,  a2,  c.b, (int)r);
        if (c.rel == Relation::GreaterEq || c.rel == Relation::Equal)
            addPlane(-a0, -a1, -a2, -c.b, (int)r);
    }
    // axis bounds and safety box
    addPlane(-1, 0, 0, 0, -1);  addPlane(0, -1, 0, 0, -1);  addPlane(0, 0, -1, 0, -1);
    addPlane( 1, 0, 0, box, -1); addPlane(0, 1, 0, box, -1); addPlane(0, 0, 1, box, -1);

    const int np = (int)P.planes.size();

    // ---- vertex enumeration: all plane triples ----------------------------
    for (int i = 0; i < np; ++i)
        for (int j = i + 1; j < np; ++j)
            for (int k = j + 1; k < np; ++k)
            {
                Pt3 v;
                if (!intersect3Planes(P.planes[i], P.planes[j], P.planes[k], v))
                    continue;

                bool inside = true;
                for (int t = 0; t < np && inside; ++t)
                    if (P.planes[t].eval(v) > tol * (1.0 + std::fabs(P.planes[t].rhs)))
                        inside = false;
                if (!inside)
                    continue;

                bool dup = false;
                for (const Pt3& w : P.verts)
                    if (norm3(v - w) < 1e-6 * (1.0 + norm3(v)))
                    {
                        dup = true;
                        break;
                    }
                if (!dup)
                    P.verts.push_back(v);
            }

    if (P.verts.size() < 3)
        return P;

    // ---- faces: incident vertices per plane, sorted by angle --------------
    for (int pl = 0; pl < np; ++pl)
    {
        const Plane3& plane = P.planes[pl];
        std::vector<int> idx;
        for (int v = 0; v < (int)P.verts.size(); ++v)
            if (std::fabs(plane.eval(P.verts[v])) < 1e-6 * (1.0 + std::fabs(plane.rhs)))
                idx.push_back(v);
        if (idx.size() < 3)
            continue;

        Pt3 cen{0,0,0};
        for (int v : idx)
            cen = cen + P.verts[v];
        cen = cen * (1.0 / (double)idx.size());

        // local 2D basis (u, w) inside the plane
        Pt3 u = std::fabs(plane.n.x) < 0.9 ? cross3(plane.n, {1,0,0})
                                           : cross3(plane.n, {0,1,0});
        u = u * (1.0 / norm3(u));
        Pt3 w = cross3(plane.n, u);

        std::vector<std::pair<double,int>> ang;
        ang.reserve(idx.size());
        for (int v : idx)
        {
            Pt3 dvec = P.verts[v] - cen;
            ang.push_back({ std::atan2(dot3(dvec, w), dot3(dvec, u)), v });
        }
        std::sort(ang.begin(), ang.end());

        Face3 f;
        f.plane    = pl;
        f.normal   = plane.n;
        f.centroid = cen;
        for (auto& a : ang)
            f.verts.push_back(a.second);
        P.faces.push_back(std::move(f));
    }
    return P;
}

} // namespace lp
