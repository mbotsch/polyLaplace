#include "PolygonSampler.h"

#include "samples.inl"

PolygonSampler::PolygonSampler(const Eigen::Matrix<double, -1, 2>& poly,
                               const double area, const double avgNumSamples)
{
    n = (int)poly.rows();

    Eigen::Matrix2d randRot;
    randRot.setRandom();
    randRot.row(0).normalize();
    randRot.row(1) =
        (randRot.row(1) - randRot.row(0) * (randRot.row(1).dot(randRot.row(0))))
            .normalized();
    randRot.setIdentity();

    double scale = sqrt(area / 4. * 1024. / avgNumSamples) /
                   sqrt(2.); ///(0.25 * 1024.) / (avgNumSamples / area);

    for (int i = 0; i < N; ++i)
    {
        Eigen::Vector2d p(scale * samples[i][0], scale * samples[i][1]);
        if (p.cwiseAbs().maxCoeff() > 1.)
            continue;
        p = randRot * p;

        // check for point in polygon
        int cnt = 0;

        for (int j = 0; j < n; ++j)
        {
            Eigen::Vector2d ev = poly.row((j + 1) % n) - poly.row(j);

            const double lambda = (p(0) - poly(j, 0)) / ev(0);
            if (lambda < 0 || lambda > 1)
                continue;

            const double mu =
                poly(j, 1) - p(1) + (ev(1) * (p(0) - poly(j, 0))) / ev(0);

            if (mu > 0)
            {
                ev.normalize();

                if (p(0) * ev(1) - p(1) * ev(0) >
                    poly(j, 0) * ev(1) - poly(j, 1) * ev(0))
                    ++cnt;
                else
                    --cnt;
            }
        }

        if (cnt)
            polySamples.push_back(p);
    }
}
