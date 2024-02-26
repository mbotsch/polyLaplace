#ifndef GLOBALENUMS_H
#define GLOBALENUMS_H

#include <string>

enum LaplaceMode2D
{
    PolySimpleLaplace = 0,
    AlexaWardetzkyLaplace = 1,
    Diamond = 2,
    deGoesLaplace = 3,
    Harmonic = 4,
};

enum InsertedPoint
{
    Centroid_ = 0,
    AreaMinimizer = 2,
    TraceMinimizer = 3,
};

enum DiffusionStep
{
    MeanEdge = 0,
    MaxEdge = 1,
    MaxDiagonal = 2
};

class SmoothingConfigs
{
public:
    explicit SmoothingConfigs(int numIters, bool fixBoundary = false,
                              bool updateQuadrics = false,
                              bool withCnum = false,
                              bool generalizedCnum = false)
        : numIters(numIters),
          fixBoundary(fixBoundary),
          updateQuadrics(updateQuadrics),
          withCnum(withCnum),
          generalizedCnum(generalizedCnum){};

    int numIters;
    bool fixBoundary;
    bool updateQuadrics;
    bool withCnum;
    bool generalizedCnum;
};

#endif
