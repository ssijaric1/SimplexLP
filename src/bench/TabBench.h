//  SimplexLP - Project #13: Revised Simplex in natID
//  Student: Šehzada Sijarić (19964)
//
//  TabBench.h
//  Benchmark tab implementing the computational comparison from the
//  proposal: a sweep of synthetic LP instances of increasing size is solved
//  by the revised simplex and the interior-point method; solve time and
//  iteration count are charted against the number of rows m, and every
//  instance can be exported to MPS for the external CLP cross-check
//  (scripts/run_clp.sh runs CLP on exactly those files).
#pragma once

#include <gui/View.h>
#include <gui/SplitterLayout.h>
#include <gui/VerticalLayout.h>
#include <gui/HorizontalLayout.h>
#include <gui/GridLayout.h>
#include <gui/GridComposer.h>
#include <gui/Label.h>
#include <gui/NumericEdit.h>
#include <gui/CheckBox.h>
#include <gui/Button.h>
#include <gui/TextEdit.h>

#include "BenchRunner.h"
#include "ViewChartCanvas.h"

class TabBench : public gui::View
{
    // ---- left: parameters ---------------------------------------------------
    class Controls : public gui::View
    {
    public:
        gui::Label _lblMStart, _lblMEnd, _lblSteps, _lblRatio, _lblNnz, _lblSeed;
        gui::NumericEdit _neMStart, _neMEnd, _neSteps, _neRatio, _neNnz, _neSeed;
        gui::CheckBox _cbSimplex, _cbIPM, _cbExport;
        gui::Button   _btnRun;
        gui::TextEdit _log;

        gui::GridLayout     _grid;
        gui::VerticalLayout _vl;

        Controls()
        : _lblMStart(tr("lblMStart"))
        , _lblMEnd(tr("lblMEnd"))
        , _lblSteps(tr("lblSteps"))
        , _lblRatio(tr("lblRatio"))
        , _lblNnz(tr("lblNnz"))
        , _lblSeed(tr("lblSeed"))
        , _neMStart(td::int4)
        , _neMEnd(td::int4)
        , _neSteps(td::int4)
        , _neRatio(td::real8)
        , _neNnz(td::int4)
        , _neSeed(td::int4)
        , _cbSimplex(tr("cbRunSimplex"))
        , _cbIPM(tr("cbRunIPM"))
        , _cbExport(tr("cbExportMPS"))
        , _btnRun(tr("btnRun"))
        , _log(gui::TextEdit::HorizontalScroll::Yes,
               gui::TextEdit::Events::DoNotSend, true)
        , _grid(6, 2)
        , _vl(6)
        {
            _neMStart.setValue(td::Variant((td::INT4)20));
            _neMEnd.setValue(td::Variant((td::INT4)200));
            _neSteps.setValue(td::Variant((td::INT4)6));
            _neRatio.setValue(td::Variant(2.5));
            _neNnz.setValue(td::Variant((td::INT4)4));
            _neSeed.setValue(td::Variant((td::INT4)42));

            _cbSimplex.setChecked(true);
            _cbIPM.setChecked(true);
            _cbExport.setChecked(false);

            gui::GridComposer gc(_grid);
            gc.appendRow(_lblMStart) << _neMStart;
            gc.appendRow(_lblMEnd)   << _neMEnd;
            gc.appendRow(_lblSteps)  << _neSteps;
            gc.appendRow(_lblRatio)  << _neRatio;
            gc.appendRow(_lblNnz)    << _neNnz;
            gc.appendRow(_lblSeed)   << _neSeed;

            _vl << _grid << _cbSimplex << _cbIPM << _cbExport << _btnRun << _log;
            setLayout(&_vl);
        }

        template <typename T>
        T num(const gui::NumericEdit& ne, T fallback) const
        {
            T v = fallback;
            ne.getValue(v);
            return v;
        }
    };

    Controls            _controls;
    ViewChartCanvas     _chartTime;
    ViewChartCanvas     _chartIters;
    gui::VerticalLayout _vlCharts;
    gui::View           _chartHost;
    gui::SplitterLayout _splitter;

    BenchRunner _runner;
    gui::thread::MainThreadSharedFunction1 _uiTick;

public:
    TabBench()
    : _vlCharts(2)
    , _chartHost(0, 0, 0, 0)
    , _splitter(gui::SplitterLayout::Orientation::Horizontal,
                gui::SplitterLayout::AuxiliaryCell::First)
    {
        setMargins(0, 0, 0, 0);

        _chartTime.setLabels(tr("chartTimeTitle"), "m", "ms");
        _chartIters.setLabels(tr("chartItersTitle"), "m", tr("chartItersY"));

        _vlCharts << _chartTime << _chartIters;
        _chartHost.setLayout(&_vlCharts);

        _splitter.setContent(_controls, _chartHost);
        setLayout(&_splitter);

        _controls._btnRun.onClick([this]() { startRun(); });

        // marshal worker progress onto the GUI thread
        _uiTick = std::make_shared<gui::thread::MainThreadFunction1>(
            [this](td::Variant v)
            {
                td::INT4 done = 0;
                v.getValue(done);
                onTick(done != 0);
            });
    }

    ~TabBench()
    {
        _runner.cancelAndJoin();
    }

protected:
    void startRun()
    {
        if (_runner.running())
            return;

        BenchParams p;
        p.mStart    = _controls.num<td::INT4>(_controls._neMStart, 20);
        p.mEnd      = _controls.num<td::INT4>(_controls._neMEnd, 200);
        p.steps     = _controls.num<td::INT4>(_controls._neSteps, 6);
        p.nOverM    = _controls.num<double>(_controls._neRatio, 2.5);
        p.nnzPerCol = _controls.num<td::INT4>(_controls._neNnz, 4);
        p.seed      = (unsigned)_controls.num<td::INT4>(_controls._neSeed, 42);
        p.runSimplex = _controls._cbSimplex.isChecked();
        p.runIPM     = _controls._cbIPM.isChecked();
        p.exportMPS  = _controls._cbExport.isChecked();

        if (p.mStart < 2)          p.mStart = 2;
        if (p.mEnd < p.mStart)     p.mEnd = p.mStart;
        if (p.steps < 1)           p.steps = 1;
        if (p.nOverM < 1.1)        p.nOverM = 1.1;
        if (p.nnzPerCol < 1)       p.nnzPerCol = 1;
        if (!p.runSimplex && !p.runIPM)
        {
            _controls._log.appendString(tr("logPickSolver"));
            _controls._log.appendString("\n");
            return;
        }

        _controls._log.setText(tr("logBenchStart"));
        _controls._log.appendString("\n");
        _controls._btnRun.disable();
        _chartTime.clearSeries();
        _chartIters.clearSeries();

        _runner.start(p, _uiTick);
    }

    void onTick(bool done)
    {
        for (const std::string& line : _runner.takeLogLines())
        {
            _controls._log.appendString(line.c_str());
            _controls._log.appendString("\n");
        }

        if (!done)
            return;

        _runner.joinIfFinished();
        _controls._btnRun.enable();
        rebuildCharts();
    }

    void rebuildCharts()
    {
        const std::vector<BenchPoint> pts = _runner.points();

        std::vector<double> xs, spxMs, ipmMs, spxIt, ipmIt;
        bool anySpx = false, anyIpm = false;
        for (const BenchPoint& p : pts)
        {
            xs.push_back((double)p.m);
            spxMs.push_back(p.spxMs);
            ipmMs.push_back(p.ipmMs);
            spxIt.push_back((double)p.spxIt);
            ipmIt.push_back((double)p.ipmIt);
            anySpx |= p.haveSpx;
            anyIpm |= p.haveIpm;
        }

        _chartTime.clearSeries();
        _chartIters.clearSeries();
        if (anySpx)
        {
            _chartTime.addSeries(xs, spxMs, td::ColorID::Crimson, tr("seriesSimplex"));
            _chartIters.addSeries(xs, spxIt, td::ColorID::Crimson, tr("seriesSimplex"));
        }
        if (anyIpm)
        {
            _chartTime.addSeries(xs, ipmMs, td::ColorID::DarkMagenta, tr("seriesIPM"));
            _chartIters.addSeries(xs, ipmIt, td::ColorID::DarkMagenta, tr("seriesIPM"));
        }
    }
};
