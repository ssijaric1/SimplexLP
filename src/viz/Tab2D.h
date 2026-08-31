//  SimplexLP - Project #13: Revised Simplex in natID
//  Student: Šehzada Sijarić (19964)
//
//  Tab2D.h
//  Two-variable visualization tab. Left: model editor + animation controls
//  (ViewVizControls). Right: View2DCanvas drawing the feasible polygon, the
//  simplex pivot path along its vertices and the interior-point trajectory.
//
//  A gui::Timer drives the "Play" animation; each tick advances the cursor
//  through the recorded pivot sequence and redraws the canvas.
#pragma once

#include <gui/View.h>
#include <gui/SplitterLayout.h>
#include <gui/Timer.h>
#include "VizModel.h"
#include "ViewVizControls.h"
#include "View2DCanvas.h"

class Tab2D : public gui::View
{
    lp::VizModel        _model;
    ViewVizControls     _controls;
    View2DCanvas        _canvas;
    gui::SplitterLayout _splitter;
    gui::Timer          _timer;

public:
    Tab2D()
    : _model(2)
    , _canvas(&_model)
    , _splitter(gui::SplitterLayout::Orientation::Horizontal,
                gui::SplitterLayout::AuxiliaryCell::First)
    , _timer(this, 0.25f, false)
    {
        setMargins(0, 0, 0, 0);
        _splitter.setContent(_controls, _canvas);
        setLayout(&_splitter);

        _controls.onSolve   = [this]() { solve(); };
        _controls.onStep    = [this]() { step(); };
        _controls.onPlay    = [this]() { togglePlay(); };
        _controls.onReset   = [this]() { reset(); };
        _controls.onShowIPM = [this](bool b) { _canvas.showIPM(b); };

        _controls.setPresets(makePresets());
    }

    ~Tab2D()
    {
        if (_timer.isRunning())
            _timer.stop();
    }

    // ----- interface used by MainWindow / toolbar --------------------------
    void solve()
    {
        stopAnimation();

        if (!_model.build(_controls.objectiveText(),
                          _controls.isMaximize(),
                          _controls.constraintsText()))
        {
            _controls.setLog(std::string(tr("errParse").c_str()) + "\n" + _model.lastError);
            _canvas.reDraw();
            return;
        }

        _model.solveBoth();
        _model.resetCursor();
        _controls.setLog(_model.summary());
        _canvas.resetView();
    }

    void step()
    {
        stopAnimation();
        if (!_model.hasPath())
        {
            solve();
            return;
        }
        if (!_model.stepCursor())
            _controls.appendLog(std::string(tr("logPathEnd").c_str()));
        _canvas.reDraw();
    }

    void togglePlay()
    {
        if (_timer.isRunning())
        {
            stopAnimation();
            return;
        }
        if (!_model.hasPath())
            solve();
        if (!_model.hasPath())
            return;

        if (_model.cursorAtEnd())
            _model.resetCursor();

        const float interval = 1.0f / (float)_controls.stepsPerSecond();
        _timer.setInterval(interval);
        _timer.start();
        _controls.setPlaying(true);
    }

    void reset()
    {
        stopAnimation();
        _model.resetCursor();
        _canvas.resetView();
    }

    bool isAnimating() const { return _timer.isRunning(); }

protected:
    void stopAnimation()
    {
        if (_timer.isRunning())
            _timer.stop();
        _controls.setPlaying(false);
    }

    bool onTimer(gui::Timer* pTimer) override
    {
        if (pTimer != &_timer)
            return false;

        if (!_model.stepCursor())
        {
            stopAnimation();
        }
        _canvas.reDraw();
        return true;
    }

    // ----- classic teaching examples ---------------------------------------
    static std::vector<ViewVizControls::Preset> makePresets()
    {
        std::vector<ViewVizControls::Preset> v;

        v.push_back({"Production (max, Dantzig)", true,
            "3x + 5y",
            "x <= 4\n2y <= 12\n3x + 2y <= 18"});

        v.push_back({"Diet (min, phase I)", false,
            "2x + 3y",
            "x + y >= 4\nx + 3y >= 6"});

        v.push_back({"Degenerate vertex", true,
            "x + y",
            "x <= 4\ny <= 4\nx + y <= 8"});

        v.push_back({"Unbounded", true,
            "x + y",
            "-x + y <= 2\nx - y <= 2"});

        v.push_back({"Infeasible", true,
            "x + y",
            "x + y <= 2\nx + y >= 5"});

        v.push_back({"Equality constraint", false,
            "x + 2y",
            "x + y = 4\nx <= 3"});

        return v;
    }
};
