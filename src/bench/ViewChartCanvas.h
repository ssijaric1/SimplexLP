//  SimplexLP - Project #13: Revised Simplex in natID
//  Student: Šehzada Sijarić (19964)
//
//  ViewChartCanvas.h
//  Small, self-contained line-chart widget on a plain gui::Canvas: axes
//  with "nice" tick steps, any number of (x, y) series with point markers,
//  and a legend. It deliberately uses only the Canvas / Shape /
//  DrawableString APIs that the SDK examples exercise, so the benchmark tab
//  carries no extra library dependency. (The SDK also ships natPlot's
//  gui::plot::View, noted as an alternative in the implementation report.)
#pragma once

#include <gui/Canvas.h>
#include <gui/Shape.h>
#include <gui/DrawableString.h>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

class ViewChartCanvas : public gui::Canvas
{
public:
    struct Series
    {
        std::vector<double> x, y;
        td::ColorID color = td::ColorID::Crimson;
        td::String  name;
    };

private:
    std::vector<Series> _series;
    td::String _title, _xName, _yName;
    gui::Size  _size {600, 260};

public:
    ViewChartCanvas() : gui::Canvas({})
    {
        enableResizeEvent(true);
    }

    void setLabels(const td::String& title, const td::String& xName,
                   const td::String& yName)
    {
        _title = title;
        _xName = xName;
        _yName = yName;
    }

    void clearSeries()
    {
        _series.clear();
        reDraw();
    }

    void addSeries(const std::vector<double>& x, const std::vector<double>& y,
                   td::ColorID color, const td::String& name)
    {
        Series s;
        s.x = x;
        s.y = y;
        s.color = color;
        s.name  = name;
        _series.push_back(std::move(s));
        reDraw();
    }

protected:
    void onResize(const gui::Size& newSize) override { _size = newSize; }

    static double niceStep(double rough)
    {
        if (rough <= 0.0 || !std::isfinite(rough))
            return 1.0;
        const double p = std::pow(10.0, std::floor(std::log10(rough)));
        const double r = rough / p;
        if (r < 1.5) return 1.0 * p;
        if (r < 3.5) return 2.0 * p;
        if (r < 7.5) return 5.0 * p;
        return 10.0 * p;
    }

    void onDraw(const gui::Rect& /*rect*/) override
    {
        const double L = 64, R = 18, T = 30, B = 42;
        const double w = (double)_size.width;
        const double h = (double)_size.height;
        const double plotW = std::max(10.0, w - L - R);
        const double plotH = std::max(10.0, h - T - B);

        // title
        gui::DrawableString::draw(_title, gui::Point(L, 6),
                                  gui::Font::ID::SystemBold, td::ColorID::SysText);

        // frame
        const gui::Point fr[4] = { {L, T}, {L + plotW, T},
                                   {L + plotW, T + plotH}, {L, T + plotH} };
        gui::Shape frame;
        frame.createPolygon(fr, 4, 1.0f);
        frame.drawWire(td::ColorID::Gray);

        // data ranges
        bool any = false;
        double xMin = 0, xMax = 1, yMin = 0, yMax = 1;
        for (const Series& s : _series)
        {
            for (size_t i = 0; i < s.x.size() && i < s.y.size(); ++i)
            {
                if (!any)
                {
                    xMin = xMax = s.x[i];
                    yMin = yMax = s.y[i];
                    any = true;
                }
                xMin = std::min(xMin, s.x[i]); xMax = std::max(xMax, s.x[i]);
                yMin = std::min(yMin, s.y[i]); yMax = std::max(yMax, s.y[i]);
            }
        }
        if (!any)
        {
            gui::DrawableString::draw(tr("chartNoData"),
                                      gui::Point(L + 12, T + 12),
                                      gui::Font::ID::SystemNormal,
                                      td::ColorID::Gray);
            return;
        }

        yMin = std::min(0.0, yMin);                  // anchor at zero
        if (xMax - xMin < 1e-12) xMax = xMin + 1.0;
        if (yMax - yMin < 1e-12) yMax = yMin + 1.0;
        yMax += (yMax - yMin) * 0.08;                // headroom

        auto X = [&](double v) { return L + (v - xMin) / (xMax - xMin) * plotW; };
        auto Y = [&](double v) { return T + plotH - (v - yMin) / (yMax - yMin) * plotH; };

        // ticks + grid
        const double sx = niceStep((xMax - xMin) / 6.0);
        const double sy = niceStep((yMax - yMin) / 5.0);
        td::String lbl;
        for (double v = std::ceil(xMin / sx) * sx; v <= xMax + 1e-9; v += sx)
        {
            const double px = X(v);
            gui::Shape::drawLine({px, T}, {px, T + plotH},
                                 td::ColorID::LightGray, 0.6f, td::LinePattern::Dot);
            lbl.format("%g", v);
            gui::DrawableString::draw(lbl, gui::Point(px - 8, T + plotH + 6),
                                      gui::Font::ID::SystemSmaller, td::ColorID::SysText);
        }
        for (double v = std::ceil(yMin / sy) * sy; v <= yMax + 1e-9; v += sy)
        {
            const double py = Y(v);
            gui::Shape::drawLine({L, py}, {L + plotW, py},
                                 td::ColorID::LightGray, 0.6f, td::LinePattern::Dot);
            lbl.format("%g", v);
            gui::DrawableString::draw(lbl, gui::Point(8, py - 7),
                                      gui::Font::ID::SystemSmaller, td::ColorID::SysText);
        }

        // axis names
        gui::DrawableString::draw(_xName, gui::Point(L + plotW * 0.5 - 12, h - 20),
                                  gui::Font::ID::SystemSmaller, td::ColorID::SysText);
        gui::DrawableString::draw(_yName, gui::Point(6, T - 18),
                                  gui::Font::ID::SystemSmaller, td::ColorID::SysText);

        // series
        for (const Series& s : _series)
        {
            const size_t nPts = std::min(s.x.size(), s.y.size());
            if (nPts == 0)
                continue;
            std::vector<gui::Point> pts;
            pts.reserve(nPts);
            for (size_t i = 0; i < nPts; ++i)
                pts.push_back({X(s.x[i]), Y(s.y[i])});

            if (nPts > 1)
            {
                gui::Shape line;
                line.createPolyLine(pts.data(), pts.size(), 2.2f);
                line.drawWire(s.color);
            }
            for (const gui::Point& p : pts)
            {
                gui::Shape dot;
                dot.createCircle(gui::Circle(p, 3.0), 1.0f);
                dot.drawFill(s.color);
            }
        }

        // legend (top-right corner of the plot area)
        double ly = T + 8;
        for (const Series& s : _series)
        {
            const double lx = L + plotW - 150;
            gui::Shape::drawLine({lx, ly + 6}, {lx + 26, ly + 6}, s.color, 3.0f);
            gui::DrawableString::draw(s.name, gui::Point(lx + 32, ly - 1),
                                      gui::Font::ID::SystemSmaller,
                                      td::ColorID::SysText);
            ly += 18;
        }
    }
};
