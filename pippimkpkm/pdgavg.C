#include <iostream>

using namespace std;

// takes two vectors of doubles (xi = values, dxi = errors) and
// returns a pair of doubles (x_avg, dx_avg)
pair<double, double> pdgavg(vector<double> xi, vector<double> dxi)
{
    int n = xi.size();
    
    // sum up weights and weighted numbers
    double sumwixi = 0, sumwi = 0;
    vector<double> wi;
    for (int i=0; i<n; ++i)
    {
        wi.push_back(1./(dxi[i]*dxi[i]));
        
        sumwixi += xi[i]*wi[i];
        sumwi   += wi[i]; 
    }
    
    // compute weighted average and error
    double x  = sumwixi/sumwi;
    double dx = 1./sqrt(sumwi);
    
    // determine chi2/ndf
    double chi2 = 0;
    for (int i=0; i<n; ++i) chi2 += wi[i]*(x - xi[i])*(x - xi[i]);
    
    // if chi2/ndf>1, scale error by sqrt(chi2/ndf)
    double S = 1.0;
    if (chi2/(n-1)>1)
    {
        S = sqrt(chi2/(n-1));
        cout <<"Applied scale factor S="<<S<<" since chi2/ndf="<<chi2/(n-1)<<endl;
    }
    
    return make_pair(sumwixi/sumwi, S/sqrt(sumwi));
}
