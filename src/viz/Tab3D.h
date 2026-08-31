//  SimplexLP - Project #13: Revised Simplex in natID
//  Student: Šehzada Sijarić (19964)
//
//  Tab3D.h
//  Three-variable visualization tab: same control panel as the 2D tab, the
//  canvas renders the feasible polytope and the pivot sequence in 3D.
#pragma once

#include <gui/View.h>
#include <gui/SplitterLayout.h>
#include <gui/Timer.h>
#include "VizModel.h"
#include "ViewVizControls.h"
#include "View3DCanvas.h"

class Tab3D : public gui::View
{
    lp::VizModel        _model;
    ViewVizControls     _controls;
    View3DCanvas        _canvas;
    gui::SplitterLayout _splitter;
    gui::Timer          _timer;

public:
    Tab3D()
    : _model(3)
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

    ~Tab3D()
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
        _canvas.reDraw();
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
            stopAnimation();
        _canvas.reDraw();
        return true;
    }

    // ----- classic 3-variable examples --------------------------------------
    static std::vector<ViewVizControls::Preset> makePresets()
    {
        std::vector<ViewVizControls::Preset> v;

        v.push_back({"Box corner (max)", true,
            "2x + 3y + z",
            "x <= 4\ny <= 3\nz <= 5\nx + y + z <= 9"});

        v.push_back({"Cut simplex (max)", true,
            "x + 2y + 3z",
            "x + y + z <= 6\nx + 2z <= 8\ny + z <= 5"});

        v.push_back({"Transport mix (min)", false,
            "4x + 3y + 5z",
            "x + y + z >= 5\n2x + y >= 4\ny + 2z >= 3"});

        v.push_back({"Equality slice", true,
            "x + y + z",
            "x + y + z = 5\nx <= 3\ny <= 3\nz <= 3"});

        v.push_back({"Unbounded wedge", true,
            "x + y + z",
            "x - y <= 2\ny - z <= 2"});

        return v;
    }
};
