//  SimplexLP - Project #13: Revised Simplex in natID
//  Student: Šehzada Sijarić (19964)
//
//  BenchRunner.h
//  Background benchmark engine for the comparison required by the proposal:
//  a sweep of synthetic LP instances of growing size is generated
//  (LPGenerator), each instance is solved by the revised simplex and by the
//  interior-point method, and solve time / iteration counts are collected.
//
//  Instances can be exported to MPS so the *same* files can be fed to
//  COIN-OR CLP (scripts/run_clp.sh) for the external comparison; a CSV with
//  all measurements is always written next to them.
//
//  The sweep runs on a std::thread; progress is marshalled back to the GUI
//  thread with gui::thread::asyncExecInMainThread, the same pattern the SDK
//  uses in its animation examples.
#pragma once

#include <gui/Thread.h>
#include "../core/LPGenerator.h"
#include "../core/RevisedSimplex.h"
#include "../core/InteriorPoint.h"
#include "../core/MPS.h"
#include "../core/Verify.h"

#include <thread>
#include <mutex>
#include <atomic>
#include <deque>
#include <vector>
#include <string>
#include <cstdlib>
#include <cmath>
#include <filesystem>
#include <fstream>

struct BenchParams
{
    int      mStart    = 20;
    int      mEnd      = 200;
    int      steps     = 6;
    double   nOverM    = 2.5;   // columns per row
    int      nnzPerCol = 4;
    unsigned seed      = 42;
    bool     runSimplex = true;
    bool     runIPM     = true;
    bool     exportMPS  = false;
    std::string outDir;         // filled by BenchRunner if empty
};

struct BenchPoint
{
    int m = 0, n = 0;
    bool haveSpx = false, haveIpm = false;

    lp::Status spxSt = lp::Status::NotSolved;
    lp::Status ipmSt = lp::Status::NotSolved;

    double spxMs = 0.0, ipmMs = 0.0;
    int    spxIt = 0,   ipmIt = 0;
    int    spxFact = 0;
    double spxObj = 0.0, ipmObj = 0.0;
};

class BenchRunner
{
    std::thread              _worker;
    mutable std::mutex       _mx;
    std::vector<BenchPoint>  _points;
    std::deque<std::string>  _pendingLog;
    std::atomic<bool>        _running {false};
    std::atomic<bool>        _cancel  {false};
    std::atomic<bool>        _joinable{false};

public:
    ~BenchRunner() { cancelAndJoin(); }

    bool running() const { return _running.load(); }

    void cancelAndJoin()
    {
        _cancel.store(true);
        joinIfFinished(true);
    }

    // Call from the main thread once the done-tick arrives (or forced).
    void joinIfFinished(bool force = false)
    {
        if (_joinable.load() || force)
        {
            if (_worker.joinable())
                _worker.join();
            _joinable.store(false);
        }
    }

    std::vector<BenchPoint> points() const
    {
        std::lock_guard<std::mutex> g(_mx);
        return _points;
    }

    std::vector<std::string> takeLogLines()
    {
        std::lock_guard<std::mutex> g(_mx);
        std::vector<std::string> out(_pendingLog.begin(), _pendingLog.end());
        _pendingLog.clear();
        return out;
    }

    static std::string defaultOutDir()
    {
        const char* home = std::getenv("HOME");
        if (!home)
            home = std::getenv("USERPROFILE"); // Windows
        std::string base = home ? home : ".";
        return base + "/SimplexLP_bench";
    }

    // uiTick runs on the GUI thread; param Variant carries 1 when finished.
    bool start(const BenchParams& params,
               const gui::thread::MainThreadSharedFunction1& uiTick)
    {
        if (_running.load())
            return false;

        joinIfFinished(true); // reap a previous run, if any

        {
            std::lock_guard<std::mutex> g(_mx);
            _points.clear();
            _pendingLog.clear();
        }
        _cancel.store(false);
        _running.store(true);
        _joinable.store(false);

        BenchParams p = params;
        if (p.outDir.empty())
            p.outDir = defaultOutDir();

        _worker = std::thread([this, p, uiTick]() { run(p, uiTick); });
        return true;
    }

private:
    void log(const std::string& line)
    {
        std::lock_guard<std::mutex> g(_mx);
        _pendingLog.push_back(line);
    }

    void post(const gui::thread::MainThreadSharedFunction1& uiTick, int done)
    {
        gui::thread::asyncExecInMainThread(uiTick, td::Variant((td::INT4)done));
    }

    void run(BenchParams p, gui::thread::MainThreadSharedFunction1 uiTick)
    {
        std::error_code ec;
        std::filesystem::create_directories(p.outDir, ec);
        if (ec)
            log("! cannot create output directory " + p.outDir + " : " + ec.message());
        else
            log("output directory: " + p.outDir);

        if (p.steps < 1)   p.steps = 1;
        if (p.mEnd < p.mStart) p.mEnd = p.mStart;

        for (int s = 0; s < p.steps && !_cancel.load(); ++s)
        {
            const double t = (p.steps == 1) ? 0.0 : (double)s / (double)(p.steps - 1);
            const int m = (int)std::lround(p.mStart + t * (p.mEnd - p.mStart));
            const int n = std::max(m + 1, (int)std::lround(m * p.nOverM));

            lp::GeneratorParams gp;
            gp.m = m;
            gp.n = n;
            gp.nnzPerCol = p.nnzPerCol;
            gp.seed = p.seed + (unsigned)s * 1000u;

            lp::StandardLP prob = lp::generateInstance(gp);

            BenchPoint pt;
            pt.m = m;
            pt.n = n;

            char buf[256];

            if (p.exportMPS)
            {
                std::snprintf(buf, sizeof(buf), "%s/bench_m%04d_n%05d.mps",
                              p.outDir.c_str(), m, n);
                std::string err;
                if (lp::writeMPS(prob, buf, &err))
                    log(std::string("wrote ") + buf);
                else
                    log("! MPS export failed: " + err);
            }

            if (p.runSimplex && !_cancel.load())
            {
                lp::SimplexOptions so;
                so.recordPath = false;
                lp::RevisedSimplex spx;
                lp::LPResult r = spx.solve(prob, so);
                pt.haveSpx = true;
                pt.spxSt   = r.status;
                pt.spxMs   = r.stats.msTotal;
                pt.spxIt   = r.stats.iterations;
                pt.spxFact = r.stats.factorizations;
                pt.spxObj  = r.objective;
            }

            if (p.runIPM && !_cancel.load())
            {
                lp::IPMOptions io;
                io.recordPath = false;
                lp::InteriorPoint ipm;
                lp::LPResult r = ipm.solve(prob, io);
                pt.haveIpm = true;
                pt.ipmSt   = r.status;
                pt.ipmMs   = r.stats.msTotal;
                pt.ipmIt   = r.stats.iterations;
                pt.ipmObj  = r.objective;
            }

            {
                std::lock_guard<std::mutex> g(_mx);
                _points.push_back(pt);
            }

            std::snprintf(buf, sizeof(buf),
                "m=%-4d n=%-5d | simplex: %-9s %7.1f ms %5d it | ipm: %-9s %7.1f ms %3d it",
                m, n,
                pt.haveSpx ? lp::toString(pt.spxSt) : "-", pt.spxMs, pt.spxIt,
                pt.haveIpm ? lp::toString(pt.ipmSt) : "-", pt.ipmMs, pt.ipmIt);
            log(buf);

            post(uiTick, 0);
        }

        writeCSV(p);

        if (_cancel.load())
            log("benchmark cancelled.");
        else
            log("benchmark finished.");

        _running.store(false);
        _joinable.store(true);
        post(uiTick, 1);
    }

    void writeCSV(const BenchParams& p)
    {
        const std::string path = p.outDir + "/results.csv";
        std::ofstream f(path);
        if (!f.is_open())
        {
            log("! cannot write " + path);
            return;
        }
        f << "m,n,spx_status,spx_ms,spx_iters,spx_factorizations,spx_obj,"
             "ipm_status,ipm_ms,ipm_iters,ipm_obj,obj_gap\n";

        std::vector<BenchPoint> pts = points();
        for (const BenchPoint& pt : pts)
        {
            const double gap = (pt.haveSpx && pt.haveIpm)
                             ? std::fabs(pt.spxObj - pt.ipmObj) /
                               (1.0 + std::fabs(pt.spxObj))
                             : 0.0;
            f << pt.m << ',' << pt.n << ','
              << (pt.haveSpx ? lp::toString(pt.spxSt) : "-") << ','
              << pt.spxMs << ',' << pt.spxIt << ',' << pt.spxFact << ','
              << pt.spxObj << ','
              << (pt.haveIpm ? lp::toString(pt.ipmSt) : "-") << ','
              << pt.ipmMs << ',' << pt.ipmIt << ',' << pt.ipmObj << ','
              << gap << '\n';
        }
        log("wrote " + path);
    }
};
