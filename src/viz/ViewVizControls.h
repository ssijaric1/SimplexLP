//  SimplexLP - Project #13: Revised Simplex in natID
//  Student: Šehzada Sijarić (19964)
//
//  ViewVizControls.h
//  Left-hand control panel shared by the 2D and the 3D visualization tabs:
//
//    * preset selector (fills the model editor with classic examples)
//    * objective direction (max/min) + objective expression
//    * multi-line constraint editor ("3x + 2y <= 18", one per line)
//    * Solve / Step / Play / Reset buttons + animation speed slider
//    * checkbox toggling the interior-point path overlay
//    * read-only log showing solver status, iterations and timing
//
//  The panel is purely a view: every action is forwarded to the owning tab
//  through std::function callbacks.
#pragma once

#include <gui/View.h>
#include <gui/Label.h>
#include <gui/LineEdit.h>
#include <gui/TextEdit.h>
#include <gui/ComboBox.h>
#include <gui/CheckBox.h>
#include <gui/Button.h>
#include <gui/Slider.h>
#include <gui/VerticalLayout.h>
#include <gui/HorizontalLayout.h>
#include <functional>
#include <string>
#include <vector>

class ViewVizControls : public gui::View
{
public:
    struct Preset
    {
        td::String  name;
        bool        maximize;
        std::string objective;
        std::string constraints;
    };

private:
    gui::Label    _lblPreset;
    gui::ComboBox _cbPreset;

    gui::Label    _lblObjective;
    gui::ComboBox _cbDirection;     // maximize / minimize
    gui::LineEdit _editObjective;

    gui::Label    _lblConstraints;
    gui::TextEdit _editConstraints;

    gui::CheckBox _cbShowIPM;

    gui::Button   _btnSolve;
    gui::Button   _btnStep;
    gui::Button   _btnPlay;
    gui::Button   _btnReset;

    gui::Label    _lblSpeed;
    gui::Slider   _slSpeed;

    gui::TextEdit _log;

    gui::HorizontalLayout _hlObjective;
    gui::HorizontalLayout _hlButtons;
    gui::HorizontalLayout _hlSpeed;
    gui::VerticalLayout   _vl;

    std::vector<Preset> _presets;
    bool _applyingPreset = false;

public:
    // callbacks wired by the owning tab
    std::function<void()>     onSolve;
    std::function<void()>     onStep;
    std::function<void()>     onPlay;
    std::function<void()>     onReset;
    std::function<void(bool)> onShowIPM;

public:
    ViewVizControls()
    : _lblPreset(tr("lblPreset"))
    , _lblObjective(tr("lblObjective"))
    , _lblConstraints(tr("lblConstraints"))
    , _editConstraints(gui::TextEdit::HorizontalScroll::No,
                       gui::TextEdit::Events::DoNotSend, false)
    , _cbShowIPM(tr("cbShowIPM"))
    , _btnSolve(tr("btnSolve"))
    , _btnStep(tr("btnStep"))
    , _btnPlay(tr("btnPlay"))
    , _btnReset(tr("btnReset"))
    , _lblSpeed(tr("lblSpeed"))
    , _log(gui::TextEdit::HorizontalScroll::No,
           gui::TextEdit::Events::DoNotSend, true)
    , _hlObjective(2)
    , _hlButtons(4)
    , _hlSpeed(2)
    , _vl(11)
    {
        _cbDirection.addItem(tr("dirMax"));
        _cbDirection.addItem(tr("dirMin"));
        _cbDirection.selectIndex(0);

        _cbPreset.onChangedSelection([this]()
        {
            applyPreset(_cbPreset.getSelectedIndex());
        });

        _cbShowIPM.setChecked(true);
        _cbShowIPM.onClick([this]()
        {
            if (onShowIPM)
                onShowIPM(_cbShowIPM.isChecked());
        });

        _btnSolve.onClick([this]() { if (onSolve) onSolve(); });
        _btnStep .onClick([this]() { if (onStep)  onStep();  });
        _btnPlay .onClick([this]() { if (onPlay)  onPlay();  });
        _btnReset.onClick([this]() { if (onReset) onReset(); });

        _slSpeed.setRange(1.0, 20.0, 19); // pivot steps per second
        _slSpeed.setValue(4.0);

        _hlObjective << _cbDirection << _editObjective;
        _hlButtons   << _btnSolve << _btnStep << _btnPlay << _btnReset;
        _hlSpeed     << _lblSpeed << _slSpeed;

        _vl << _lblPreset << _cbPreset;
        _vl << _lblObjective << _hlObjective;
        _vl << _lblConstraints << _editConstraints;
        _vl.appendSpace(4);
        _vl << _cbShowIPM << _hlButtons << _hlSpeed << _log;

        setLayout(&_vl);
    }

    // -------------------------------------------------------------- presets
    void setPresets(const std::vector<Preset>& presets)
    {
        _presets = presets;
        for (const Preset& p : _presets)
            _cbPreset.addItem(p.name);
        if (!_presets.empty())
        {
            _cbPreset.selectIndex(0);
            applyPreset(0);
        }
    }

    void applyPreset(int i)
    {
        if (i < 0 || i >= (int)_presets.size())
            return;
        _applyingPreset = true;
        const Preset& p = _presets[(size_t)i];
        _cbDirection.selectIndex(p.maximize ? 0 : 1);
        _editObjective.setText(p.objective.c_str());
        _editConstraints.setText(td::String(p.constraints.c_str()));
        _applyingPreset = false;
    }

    // -------------------------------------------------------------- getters
    bool isMaximize() const { return _cbDirection.getSelectedIndex() == 0; }

    std::string objectiveText() const
    {
        td::String s = _editObjective.getText();
        return std::string(s.c_str() ? s.c_str() : "");
    }

    std::string constraintsText() const
    {
        td::String s = _editConstraints.getText();
        return std::string(s.c_str() ? s.c_str() : "");
    }

    bool   showIPM() const { return _cbShowIPM.isChecked(); }
    double stepsPerSecond() const { return _slSpeed.getValue(); }

    // ----------------------------------------------------------------- log
    void setLog(const std::string& text)
    {
        _log.setText(td::String(text.c_str()));
    }

    void appendLog(const std::string& line)
    {
        std::string s = line;
        s += "\n";
        _log.appendString(s.c_str());
    }

    void setPlaying(bool playing)
    {
        _btnPlay.setTitle(playing ? tr("btnPause") : tr("btnPlay"));
    }
};
