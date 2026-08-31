//  SimplexLP - Project #13: Revised Simplex in natID
//  Student: Šehzada Sijarić (19964)
//
//  ViewSettings.h - content of the settings dialog: shows the active UI
//  language, lets the user pick another one (EN / BA translations live in
//  res/tr/), and toggles toolbar labels. Follows the SDK example pattern;
//  a language change asks for an application restart.
#pragma once
#include <gui/View.h>
#include <gui/Label.h>
#include <gui/ComboBox.h>
#include <gui/CheckBox.h>
#include <gui/LineEdit.h>
#include <gui/GridLayout.h>
#include <gui/GridComposer.h>
#include <gui/ToolBar.h>
#include <gui/Application.h>

class ViewSettings : public gui::View
{
protected:
    gui::Label    _lblLangNow;
    gui::LineEdit _leLang;
    gui::Label    _lblLangNew;
    gui::ComboBox _cmbLangs;
    gui::CheckBox _chbToolbarIconsAndLabels;
    gui::GridLayout _gl;
    gui::ToolBar* _pMainTB = nullptr;
    int _initialLangSelection;

public:
    ViewSettings()
    : _lblLangNow(tr("lblLang"))
    , _lblLangNew(tr("lblLang2"))
    , _chbToolbarIconsAndLabels(tr("chbTBIcsAndLbls"))
    , _gl(3, 2)
    {
        gui::Application* pApp = getApplication();
        auto appProperties = pApp->getProperties();
        assert(appProperties);
        td::String strTr = appProperties->getValue("translation", "EN");

        _leLang.setAsReadOnly();

        int newLangIndex = 0;
        auto& langs = getSupportedLanguages();
        auto currTranslationIndex = getTranslationLanguageIndex();
        auto& strCurrentLanguage = langs[currTranslationIndex].getDescription();
        _leLang.setText(strCurrentLanguage);

        int i = 0;
        for (const auto& lang : langs)
        {
            if (lang.getExtension() == strTr)
                newLangIndex = i;
            _cmbLangs.addItem(lang.getDescription());
            ++i;
        }

        bool showLabels = appProperties->getTBLabelVisibility(
                              mu::IAppProperties::ToolBarType::Main, true);
        _chbToolbarIconsAndLabels.setChecked(showLabels);

        _cmbLangs.selectIndex(newLangIndex);
        _initialLangSelection = newLangIndex;

        gui::GridComposer gc(_gl);
        gc.appendRow(_lblLangNow) << _leLang;
        gc.appendRow(_lblLangNew) << _cmbLangs;
        gc.appendRow(_chbToolbarIconsAndLabels, 0);
        setLayout(&_gl);

        _chbToolbarIconsAndLabels.onClick([this]()
        {
            if (_pMainTB)
                _pMainTB->showLabels(_chbToolbarIconsAndLabels.isChecked());
        });
    }

    td::String getTranslationExt()
    {
        td::String strExt;
        int currSelection = _cmbLangs.getSelectedIndex();
        if (currSelection >= 0)
        {
            auto& langs = getSupportedLanguages();
            strExt = langs[currSelection].getExtension();
        }
        return strExt;
    }

    void setMainTB(gui::ToolBar* pTB) { _pMainTB = pTB; }

    bool isRestartRequired() const
    {
        return _initialLangSelection != _cmbLangs.getSelectedIndex();
    }
};
