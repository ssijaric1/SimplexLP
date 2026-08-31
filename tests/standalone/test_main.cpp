//  SimplexLP - standalone correctness tests
//  Student: Šehzada Sijarić (19964)
//
//  Compiles the *unmodified* LP core against the shim headers (dense-LU
//  stand-in for natID's sparse solver) and verifies the mathematics:
//
//    1. textbook LPs with known optima (crash basis, phase I, max/min)
//    2. infeasible / unbounded detection
//    3. Beale's cycling example (Bland's rule safeguard)
//    4. redundant equality rows (artificials driven out / parked)
//    5. simplex vs interior-point agreement + KKT on random instances
//    6. MPS write -> read round trip
//    7. 2D polygon clipping and 3D polytope enumeration
//
//  Build:  ./build_and_run.sh   (plain g++, no SDK required)
#include <cstdio>
#include <cmath>
#include <string>

#include "../../src/core/LPCommon.h"
#include "../../src/core/LPProblem.h"
#include "../../src/core/RevisedSimplex.h"
#include "../../src/core/InteriorPoint.h"
#include "../../src/core/LPGenerator.h"
#include "../../src/core/MPS.h"
#include "../../src/core/Verify.h"
#include "../../src/core/Polytope.h"

static int gFailures = 0;
static int gChecks   = 0;

#define CHECK(cond, msg) do {                                              \
    ++gChecks;                                                             \
    if (!(cond)) {                                                         \
        ++gFailures;                                                       \
        std::printf("  FAIL %s:%d  %s\n", __FILE__, __LINE__, msg);        \
    }                                                                      \
} while (0)

static bool near(double a, double b, double tol = 1e-6)
{
    return std::fabs(a - b) <= tol * (1.0 + std::fabs(a) + std::fabs(b));
}

// ---------------------------------------------------------------------------
static void testTextbookMax()
{
    std::printf("[1] textbook max (Dantzig): max 3x+5y, x<=4, 2y<=12, 3x+2y<=18\n");
    lp::UserModel um;
    um.maximize = true;
    um.c        = {3, 5};
    um.varNames = {"x", "y"};
    um.rows.push_back({{1, 0}, lp::Relation::LessEq, 4});
    um.rows.push_back({{0, 2}, lp::Relation::LessEq, 12});
    um.rows.push_back({{3, 2}, lp::Relation::LessEq, 18});

    lp::StandardLP P = lp::toStandardForm(um);

    lp::SimplexOptions so;
    so.recordPath = true;
    lp::RevisedSimplex spx;
    lp::LPResult r = spx.solve(P, so);

    CHECK(r.status == lp::Status::Optimal, "status optimal");
    CHECK(r.stats.phase1Iters == 0, "crash basis: no phase I");
    CHECK(near(P.userObjective(&r.x[0]), 36.0), "objective = 36");
    CHECK(near(r.x[0], 2.0) && near(r.x[1], 6.0), "x=(2,6)");
    CHECK(!r.path.empty(), "path recorded");
    CHECK(near(r.path.front().x[0], 0.0) && near(r.path.front().x[1], 0.0),
          "path starts at origin (slack basis)");

    lp::KKTReport kkt = lp::verifyKKT(P, r.x, r.y);
    CHECK(kkt.ok(1e-6), "KKT satisfied");

    // phase-II objective must be non-increasing (standard min form)
    for (size_t i = 1; i < r.path.size(); ++i)
        CHECK(r.path[i].objective <= r.path[i - 1].objective + 1e-9,
              "monotone objective");
}

// ---------------------------------------------------------------------------
static void testPhase1Min()
{
    std::printf("[2] phase I needed: min 2x+3y, x+y>=10, x<=8, y<=8\n");
    lp::UserModel um;
    um.maximize = false;
    um.c        = {2, 3};
    um.rows.push_back({{1, 1}, lp::Relation::GreaterEq, 10});
    um.rows.push_back({{1, 0}, lp::Relation::LessEq, 8});
    um.rows.push_back({{0, 1}, lp::Relation::LessEq, 8});

    lp::StandardLP P = lp::toStandardForm(um);
    lp::RevisedSimplex spx;
    lp::LPResult r = spx.solve(P);

    CHECK(r.status == lp::Status::Optimal, "status optimal");
    CHECK(r.stats.phase1Iters > 0, "phase I ran");
    CHECK(near(P.userObjective(&r.x[0]), 22.0), "objective = 22");
    CHECK(near(r.x[0], 8.0) && near(r.x[1], 2.0), "x=(8,2)");

    lp::KKTReport kkt = lp::verifyKKT(P, r.x, r.y);
    CHECK(kkt.ok(1e-6), "KKT satisfied");
}

// ---------------------------------------------------------------------------
static void testEquality()
{
    std::printf("[3] equality rows: min x+2y, x+y=4, x<=3\n");
    lp::UserModel um;
    um.c = {1, 2};
    um.rows.push_back({{1, 1}, lp::Relation::Equal, 4});
    um.rows.push_back({{1, 0}, lp::Relation::LessEq, 3});

    lp::StandardLP P = lp::toStandardForm(um);
    lp::RevisedSimplex spx;
    lp::LPResult r = spx.solve(P);

    CHECK(r.status == lp::Status::Optimal, "status optimal");
    CHECK(near(P.userObjective(&r.x[0]), 5.0), "objective = 5 at (3,1)");
    CHECK(near(r.x[0], 3.0) && near(r.x[1], 1.0), "x=(3,1)");
}

// ---------------------------------------------------------------------------
static void testInfeasible()
{
    std::printf("[4] infeasible: x<=1 and x>=3\n");
    lp::UserModel um;
    um.c = {1};
    um.rows.push_back({{1}, lp::Relation::LessEq, 1});
    um.rows.push_back({{1}, lp::Relation::GreaterEq, 3});

    lp::StandardLP P = lp::toStandardForm(um);
    lp::RevisedSimplex spx;
    lp::LPResult r = spx.solve(P);
    CHECK(r.status == lp::Status::Infeasible, "detected infeasible");
}

// ---------------------------------------------------------------------------
static void testUnbounded()
{
    std::printf("[5] unbounded: max x+y, x-y<=1\n");
    lp::UserModel um;
    um.maximize = true;
    um.c = {1, 1};
    um.rows.push_back({{1, -1}, lp::Relation::LessEq, 1});

    lp::StandardLP P = lp::toStandardForm(um);
    lp::RevisedSimplex spx;
    lp::LPResult r = spx.solve(P);
    CHECK(r.status == lp::Status::Unbounded, "detected unbounded");
}

// ---------------------------------------------------------------------------
static void testBealeCycling()
{
    std::printf("[6] Beale cycling example terminates (Bland safeguard)\n");
    // min -3/4 x1 + 150 x2 - 1/50 x3 + 6 x4
    //     1/4 x1 -  60 x2 - 1/25 x3 + 9 x4 <= 0
    //     1/2 x1 -  90 x2 - 1/50 x3 + 3 x4 <= 0
    //                              x3      <= 1
    // optimum -1/20 = -0.05
    lp::UserModel um;
    um.c = {-0.75, 150.0, -0.02, 6.0};
    um.rows.push_back({{0.25,  -60.0, -1.0/25.0, 9.0}, lp::Relation::LessEq, 0});
    um.rows.push_back({{0.5,   -90.0, -1.0/50.0, 3.0}, lp::Relation::LessEq, 0});
    um.rows.push_back({{0.0,     0.0,  1.0,      0.0}, lp::Relation::LessEq, 1});

    lp::StandardLP P = lp::toStandardForm(um);
    lp::SimplexOptions so;
    so.blandAfter = 6; // make the safeguard kick in early
    lp::RevisedSimplex spx;
    lp::LPResult r = spx.solve(P, so);

    CHECK(r.status == lp::Status::Optimal, "status optimal (no cycling)");
    CHECK(near(P.userObjective(&r.x[0]), -0.05, 1e-7), "objective = -0.05");
}

// ---------------------------------------------------------------------------
static void testRedundantRows()
{
    std::printf("[7] redundant duplicated equality rows\n");
    lp::UserModel um;
    um.c = {1, 1};
    um.rows.push_back({{1, 1}, lp::Relation::Equal, 2});
    um.rows.push_back({{1, 1}, lp::Relation::Equal, 2}); // duplicate
    um.rows.push_back({{2, 2}, lp::Relation::Equal, 4}); // scaled duplicate

    lp::StandardLP P = lp::toStandardForm(um);
    lp::RevisedSimplex spx;
    lp::LPResult r = spx.solve(P);

    CHECK(r.status == lp::Status::Optimal, "status optimal");
    CHECK(near(P.userObjective(&r.x[0]), 2.0), "objective = 2");
    CHECK(near(r.x[0] + r.x[1], 2.0), "constraint satisfied");
}

// ---------------------------------------------------------------------------
static void testSimplexVsIPM()
{
    std::printf("[8] simplex vs interior point on synthetic instances\n");
    struct Cfg { int m, n, nnz; unsigned seed; };
    const Cfg cfgs[] = {
        {  8,  20, 3, 1 }, {  8,  20, 3, 2 }, {  8,  20, 3, 3 },
        { 20,  50, 4, 7 }, { 20,  50, 4, 8 },
        { 40, 100, 5, 11 }, { 40, 100, 5, 12 },
        { 60, 150, 6, 21 },
    };

    for (const Cfg& cf : cfgs)
    {
        lp::GeneratorParams gp;
        gp.m = cf.m; gp.n = cf.n; gp.nnzPerCol = cf.nnz; gp.seed = cf.seed;
        lp::StandardLP P = lp::generateInstance(gp);

        lp::RevisedSimplex spx;
        lp::LPResult rs = spx.solve(P);

        lp::InteriorPoint ipm;
        lp::IPMOptions io;
        lp::LPResult ri = ipm.solve(P, io);

        char tag[96];
        std::snprintf(tag, sizeof(tag), "m=%d n=%d seed=%u", cf.m, cf.n, cf.seed);

        CHECK(rs.status == lp::Status::Optimal,
              (std::string("simplex optimal  ") + tag).c_str());
        CHECK(ri.status == lp::Status::Optimal,
              (std::string("ipm optimal      ") + tag).c_str());

        if (rs.status == lp::Status::Optimal && ri.status == lp::Status::Optimal)
        {
            const double rel = std::fabs(rs.objective - ri.objective) /
                               (1.0 + std::fabs(rs.objective));
            CHECK(rel < 1e-5,
                  (std::string("objectives agree ") + tag).c_str());

            lp::KKTReport ks = lp::verifyKKT(P, rs.x, rs.y);
            lp::KKTReport ki = lp::verifyKKT(P, ri.x, ri.y);
            CHECK(ks.ok(1e-5), (std::string("simplex KKT ") + tag).c_str());
            CHECK(ki.ok(1e-4), (std::string("ipm KKT     ") + tag).c_str());

            std::printf("    %-22s  obj=%.8g  spxIt=%-4d ipmIt=%-3d gap=%.1e\n",
                        tag, rs.objective, rs.stats.iterations,
                        ri.stats.iterations, rel);
        }
    }
}

// ---------------------------------------------------------------------------
static void testMPSRoundTrip()
{
    std::printf("[9] MPS write -> read round trip\n");
    lp::GeneratorParams gp;
    gp.m = 15; gp.n = 40; gp.nnzPerCol = 4; gp.seed = 99;
    lp::StandardLP P1 = lp::generateInstance(gp);

    std::string err;
    CHECK(lp::writeMPS(P1, "/tmp/simplexlp_test.mps", &err), err.c_str());

    lp::StandardLP P2;
    CHECK(lp::readMPS("/tmp/simplexlp_test.mps", P2, &err), err.c_str());

    CHECK(P2.rows() == P1.rows() && P2.colsTotal() == P1.colsTotal(),
          "dimensions preserved");
    CHECK(P2.A.nnz() == P1.A.nnz(), "nnz preserved");

    lp::RevisedSimplex spx;
    lp::LPResult r1 = spx.solve(P1);
    lp::LPResult r2 = spx.solve(P2);
    CHECK(r1.status == lp::Status::Optimal && r2.status == lp::Status::Optimal,
          "both solvable");
    CHECK(near(r1.objective, r2.objective, 1e-9), "objectives identical");
}

// ---------------------------------------------------------------------------
static void testPolygon2D()
{
    std::printf("[10] 2D feasible polygon (Dantzig example)\n");
    lp::UserModel um;
    um.maximize = true;
    um.c = {3, 5};
    um.rows.push_back({{1, 0}, lp::Relation::LessEq, 4});
    um.rows.push_back({{0, 2}, lp::Relation::LessEq, 12});
    um.rows.push_back({{3, 2}, lp::Relation::LessEq, 18});

    std::vector<lp::Pt2> poly = lp::feasiblePolygon2D(um, 100.0);
    CHECK(poly.size() == 5, "5 corner points");

    const lp::Pt2 expect[5] = {{0,0},{4,0},{4,3},{2,6},{0,6}};
    for (const lp::Pt2& e : expect)
    {
        bool found = false;
        for (const lp::Pt2& p : poly)
            if (near(p.x, e.x, 1e-9) && near(p.y, e.y, 1e-9))
                found = true;
        CHECK(found, "expected vertex present");
    }

    lp::Pt2 isect;
    CHECK(lp::intersect2D(1, 0, 4, 3, 2, 18, isect), "2x2 dense solve");
    CHECK(near(isect.x, 4.0) && near(isect.y, 3.0), "intersection (4,3)");
}

// ---------------------------------------------------------------------------
static void testPolytope3D()
{
    std::printf("[11] 3D polytope (unit cube + simplex corner)\n");
    {
        lp::UserModel um;
        um.c = {1, 1, 1};
        um.rows.push_back({{1,0,0}, lp::Relation::LessEq, 1});
        um.rows.push_back({{0,1,0}, lp::Relation::LessEq, 1});
        um.rows.push_back({{0,0,1}, lp::Relation::LessEq, 1});

        lp::Polytope3 P = lp::buildPolytope3D(um, 10.0);
        CHECK(P.verts.size() == 8, "cube: 8 vertices");
        CHECK(P.faces.size() == 6, "cube: 6 faces");
    }
    {
        lp::UserModel um; // x+y+z <= 1 corner simplex
        um.c = {1, 1, 1};
        um.rows.push_back({{1,1,1}, lp::Relation::LessEq, 1});

        lp::Polytope3 P = lp::buildPolytope3D(um, 10.0);
        CHECK(P.verts.size() == 4, "tetrahedron: 4 vertices");
        CHECK(P.faces.size() == 4, "tetrahedron: 4 faces");
    }
}

// ---------------------------------------------------------------------------
static void testIPMPath()
{
    std::printf("[12] IPM path recording on the 2D demo\n");
    lp::UserModel um;
    um.maximize = true;
    um.c = {3, 5};
    um.rows.push_back({{1, 0}, lp::Relation::LessEq, 4});
    um.rows.push_back({{0, 2}, lp::Relation::LessEq, 12});
    um.rows.push_back({{3, 2}, lp::Relation::LessEq, 18});

    lp::StandardLP P = lp::toStandardForm(um);
    lp::IPMOptions io;
    io.recordPath = true;
    lp::InteriorPoint ipm;
    lp::LPResult r = ipm.solve(P, io);

    CHECK(r.status == lp::Status::Optimal, "ipm optimal");
    CHECK(near(P.userObjective(&r.x[0]), 36.0, 1e-4), "ipm objective = 36");
    CHECK(r.path.size() >= 3, "ipm path recorded");
    // last iterate close to the simplex optimum (2,6)
    const lp::IterSnapshot& last = r.path.back();
    CHECK(near(last.x[0], 2.0, 1e-3) && near(last.x[1], 6.0, 1e-3),
          "ipm converges to (2,6)");
}

// ---------------------------------------------------------------------------
int main()
{
    std::printf("SimplexLP core test suite (standalone shims)\n");
    std::printf("=============================================\n");

    testTextbookMax();
    testPhase1Min();
    testEquality();
    testInfeasible();
    testUnbounded();
    testBealeCycling();
    testRedundantRows();
    testSimplexVsIPM();
    testMPSRoundTrip();
    testPolygon2D();
    testPolytope3D();
    testIPMPath();

    std::printf("=============================================\n");
    std::printf("%d checks, %d failures\n", gChecks, gFailures);
    return gFailures == 0 ? 0 : 1;
}
