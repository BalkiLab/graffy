#include "statistics.h"

void CDLib::normalize_data(vector<double>& samples) {
    if (samples.size() == 0)
        return;
    if (samples.size() == 1) {
        samples[0] = 0;
        return;
    }
    double sum = accumulate(samples.begin(), samples.end(), 0.0);
    for (id_type i = 0; i < samples.size(); i++) {
        samples[i] /= sum;
    }
}

bool CDLib::is_distribution(const vector<double>& distribution) {
    double sum = 0;
    vector<double>::const_iterator it;
    for (it = distribution.begin(); it != distribution.end(); it++)
        sum += *it;
    if ((sum >= 0.9) && (sum <= 1.1)) // Allowing some error.
        return true;
    else
        return false;
}

double CDLib::distribution_entropy(const vector<double>& distribution) {
    double entropy = 0;
    if (is_distribution(distribution)) {
        for (unsigned int i = 0; i < distribution.size(); i++) {
            if (distribution[i] != 0) {
                entropy += distribution[i] * (log(distribution[i]) / log(2));
            }
        }
        entropy *= -1;
    }
    return entropy;
}

double CDLib::kl_divergence(const vector<double>& distribution1, const vector<double>& distribution2) {
    //        KL-Divergence is not a symmetric fucnction.
    if (!(is_distribution(distribution1) && is_distribution(distribution2)))
        return -9999; // Reporting invalid distribution error
    double sum = 0;
    id_type min = (distribution1.size() <= distribution2.size()) ? distribution1.size() : distribution2.size();
    for (id_type i = 0; i < min; i++)
        if ((distribution1[i] != 0) && (distribution2[i] != 0))
            sum += log(distribution1[i] / distribution2[i]) * distribution1[i];
    return sum;
}

double CDLib::kl_divergence_symmetric(const vector<double>& distribution1, const vector<double>& distribution2) {
    return (kl_divergence(distribution1, distribution2) + kl_divergence(distribution2, distribution1));
}

double CDLib::bhattacharyya_distance(const vector<double>& distribution1, const vector<double>& distribution2) {
    if (!(is_distribution(distribution1) && is_distribution(distribution2)))
        return -9999; // Reporting invalid distribution error
    double sum = 0;
    id_type min = (distribution1.size() <= distribution2.size()) ? distribution1.size() : distribution2.size();
    for (id_type i = 0; i < min; i++)
        sum += sqrt(distribution1[i] * distribution2[i]);
    return sum;
}

double CDLib::hellinger_distance(const vector<double>& distribution1, const vector<double>& distribution2) {
    return sqrt(1 - bhattacharyya_distance(distribution1, distribution2));
}

