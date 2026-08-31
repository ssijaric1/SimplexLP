//  SimplexLP - Project #13: Revised Simplex in natID
//  Student: Šehzada Sijarić (19964)
//
//  LinearExprParser.h
//  Tiny parser for the constraint/objective text boxes of the 2D/3D tabs.
//
//  Accepted syntax (whitespace free-form, '*' optional):
//      term      :=  [+|-] [number] [ '*' ] [variable]
//      variable  :=  x | y | z | x1 | x2 | x3   (case-insensitive)
//      expr      :=  term { (+|-) term }
//      line      :=  expr (<=|>=|=|<|>) number-expr
//
//  Examples:   3x + 2y <= 18      -x1+4 x2 >= 7      x + y + z = 10
//  Pure geometry/string code - no gui dependencies, covered by tests.
#pragma once

#include "../core/LPProblem.h"
#include <cctype>
#include <sstream>

namespace lp
{

class LinearExprParser
{
public:
    // dim = 2 or 3 (number of structural variables)
    explicit LinearExprParser(int dim) : _dim(dim) {}

    const std::string& error() const { return _err; }

    // "3x + 2y" -> coefficient vector; returns false + error() on bad input
    bool parseExpr(const std::string& text, std::vector<double>& coef, double& constant)
    {
        coef.assign((size_t)_dim, 0.0);
        constant = 0.0;
        _err.clear();

        const char* p   = text.c_str();
        bool anyTerm    = false;

        skipWs(p);
        while (*p)
        {
            double sign = 1.0;
            while (*p == '+' || *p == '-')
            {
                if (*p == '-')
                    sign = -sign;
                ++p;
                skipWs(p);
            }
            if (!*p)
            {
                _err = "dangling sign at end of expression";
                return false;
            }

            double value   = 1.0;
            bool   hasNum  = false;
            if (std::isdigit((unsigned char)*p) || *p == '.')
            {
                char* end = nullptr;
                value = std::strtod(p, &end);
                if (end == p)
                {
                    _err = "cannot read number near '" + std::string(p).substr(0, 8) + "'";
                    return false;
                }
                p = end;
                hasNum = true;
                skipWs(p);
                if (*p == '*')
                {
                    ++p;
                    skipWs(p);
                }
            }

            int var = parseVarName(p);
            if (var >= 0)
            {
                if (var >= _dim)
                {
                    _err = "variable index out of range for this tab";
                    return false;
                }
                coef[(size_t)var] += sign * value;
            }
            else if (hasNum)
            {
                constant += sign * value;
            }
            else
            {
                _err = "unexpected character '" + std::string(1, *p) + "'";
                return false;
            }
            anyTerm = true;
            skipWs(p);
        }

        if (!anyTerm)
        {
            _err = "empty expression";
            return false;
        }
        return true;
    }

    // full constraint line "expr REL expr"
    bool parseConstraint(const std::string& line, IneqConstraint& out)
    {
        _err.clear();
        size_t pos = std::string::npos;
        Relation rel = Relation::LessEq;
        size_t opLen = 2;

        if      ((pos = line.find("<=")) != std::string::npos) rel = Relation::LessEq;
        else if ((pos = line.find(">=")) != std::string::npos) rel = Relation::GreaterEq;
        else if ((pos = line.find("=<")) != std::string::npos) rel = Relation::LessEq;
        else if ((pos = line.find("=>")) != std::string::npos) rel = Relation::GreaterEq;
        else if ((pos = line.find('<'))  != std::string::npos) { rel = Relation::LessEq;    opLen = 1; }
        else if ((pos = line.find('>'))  != std::string::npos) { rel = Relation::GreaterEq; opLen = 1; }
        else if ((pos = line.find('='))  != std::string::npos) { rel = Relation::Equal;     opLen = 1; }
        else
        {
            _err = "constraint needs a relation (<=, >= or =)";
            return false;
        }

        std::vector<double> lhs, rhsCoef;
        double lhsConst = 0.0, rhsConst = 0.0;
        if (!parseExpr(line.substr(0, pos), lhs, lhsConst))
            return false;
        if (!parseExpr(line.substr(pos + opLen), rhsCoef, rhsConst))
            return false;

        // move variables left, constants right:  lhs - rhs REL rhsConst - lhsConst
        out.a.assign((size_t)_dim, 0.0);
        bool anyVar = false;
        for (int j = 0; j < _dim; ++j)
        {
            out.a[(size_t)j] = lhs[(size_t)j] - rhsCoef[(size_t)j];
            if (out.a[(size_t)j] != 0.0)
                anyVar = true;
        }
        out.rel = rel;
        out.b   = rhsConst - lhsConst;

        if (!anyVar)
        {
            _err = "constraint contains no variable";
            return false;
        }
        return true;
    }

    // Multi-line constraint block; empty lines and '#' comments are skipped.
    bool parseModel(const std::string& objective, bool maximize,
                    const std::string& constraintBlock, UserModel& out)
    {
        out = UserModel();
        out.maximize = maximize;
        out.varNames = (_dim == 2) ? std::vector<std::string>{"x", "y"}
                                   : std::vector<std::string>{"x", "y", "z"};

        double objConst = 0.0;
        if (!parseExpr(objective, out.c, objConst))
        {
            _err = "objective: " + _err;
            return false;
        }

        std::istringstream is(constraintBlock);
        std::string line;
        int lineNo = 0;
        while (std::getline(is, line))
        {
            ++lineNo;
            // trim + skip comments / blanks
            size_t h = line.find('#');
            if (h != std::string::npos)
                line = line.substr(0, h);
            bool blank = true;
            for (char ch : line)
                if (!std::isspace((unsigned char)ch))
                    blank = false;
            if (blank)
                continue;

            IneqConstraint c;
            if (!parseConstraint(line, c))
            {
                _err = "line " + std::to_string(lineNo) + ": " + _err;
                return false;
            }
            out.rows.push_back(std::move(c));
        }

        if (out.rows.empty())
        {
            _err = "no constraints entered";
            return false;
        }
        return true;
    }

private:
    int _dim;
    std::string _err;

    static void skipWs(const char*& p)
    {
        while (*p == ' ' || *p == '\t' || *p == '\r')
            ++p;
    }

    // returns variable index or -1; advances p on success
    int parseVarName(const char*& p)
    {
        const char c = (char)std::tolower((unsigned char)*p);
        if (c != 'x' && c != 'y' && c != 'z')
            return -1;
        ++p;
        if (c == 'y') return 1;
        if (c == 'z') return 2;
        // 'x' optionally followed by an index digit: x1, x2, x3
        if (*p >= '1' && *p <= '3')
        {
            const int idx = *p - '1';
            ++p;
            return idx;
        }
        return 0;
    }
};

} // namespace lp
