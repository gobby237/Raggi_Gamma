#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <iostream>
#include <boost/math/distributions/chi_squared.hpp>
using namespace std;
double normal_cdf(double z) {
    return 0.5 * std::erfc(-z / std::sqrt(2));
}
// Calcolo p-value per test z a due code
double p_value_gaussiano(double x1, double mu, double sigma) {
    double z = (x1 - mu) / sigma;
    return 2 * (1 - normal_cdf(std::fabs(z)));
}
int main() {

    string nome; 

    cout << "Inserisci nome dati: "; 
    cin >> nome; 

    ofstream out("p_value.txt", ios::app); 

    ifstream in(nome);
    if (!in.is_open()) {
        cerr << "Errore: impossibile aprire il file test.txt" << endl;
        return 1;
    }



    int dof = 30; 
    boost::math::chi_squared_distribution<double> dist(dof);

    
    vector<double> o; 

    vector<double> x;
    vector<double> y;

    double value, val;
    double min1, max1, min2, max2;

    

    // cout << "Inserire intervalli (min1 max1 min2 max2): ";
    // cin >> min1 >> max1 >> min2 >> max2;

    min1 = 320; 
    max1 = 335; 
    min2 = 390; 
    max2 = 405; 

    // Lettura e selezione dati
    while (in >> value >> val) {
    	o.push_back(val); 
        if ((min1 <= value && value <= max1) || (min2 <= value && value <= max2)) {
            x.push_back(value - 0.5);
            y.push_back(val);
        }
    }

    if (x.empty()) {
        cerr << "Errore: nessun dato negli intervalli specificati" << endl;
        return 1;
    }

    // Calcolo regressione ponderata (s_i^2 = y_i)
    double sum_w = 0.0;      // Sum(1/y_i)
    double sum_wx = 0.0;      // Sum(x_i/y_i)
    double sum_wy = 0.0;      // Sum(y_i/y_i) = Sum(1)
    double sum_wx2 = 0.0;     // Sum(x_i^2/y_i)
    double sum_wxy = 0.0;     // Sum(x_i*y_i/y_i) = Sum(x_i)
    double chiq = 0; 

    for (size_t i = 0; i < x.size(); i++) {
        double weight = 1.0 / y[i];
        sum_w += weight;
        sum_wx += x[i] * weight;
        sum_wy += 1.0;        // y[i]/y[i] = 1
        sum_wx2 += x[i] * x[i] * weight;
        sum_wxy += x[i];      // x[i]*y[i]/y[i] = x[i]
        
        
        
    }

    double delta = sum_w * sum_wx2 - sum_wx * sum_wx;

    // Calcolo parametri corretti
    double a = (sum_wx2 * sum_wy - sum_wx * sum_wxy) / delta;
    double b = (sum_w * sum_wxy - sum_wx * sum_wy) / delta;

    // Calcolo incertezze
    double sigma_a = sqrt(sum_wx2 / delta);
    double sigma_b = sqrt(sum_w / delta);
    
    double var_a = sum_wx2 / delta; 
    double var_b = sum_w / delta; 

    cout << "a = " << a << " pm " << sigma_a << endl;
    cout << "b = " << b << " pm " << sigma_b << endl;
    
    double F = 0; 
    double deltax = 385-340; 
    double mediax = (385 + 340)/2; 
    F = deltax*(a+(b*mediax));
	
	cout << "area sotto la retta tra 340 e 385: " << F << endl; 
	
	double Y = 0; 
	
	for  (int i = 340; i < 386; i ++)
	{
		Y += o.at(i); 
	}
	
	cout << "Area Y: " << Y << endl; 
	
	
	cout << "Area per differenza: " << Y-F << endl; 

	
	double S = Y - F; 
	
	double incertezza_S = 0;
	
	double var_Y = Y; 
		
	double var_F = pow(deltax,2)*(var_a + pow(mediax,2)*(var_b) - mediax*(2*sum_wx/delta)); 
	 
	
	double lambda = (S)/(sqrt(var_F + var_Y)); 

    cout << "Incertezza di S: " << sqrt(var_F + var_Y) << endl; 
	
	cout << "Lambda (compatibilità della gaussiana): " << lambda << endl; 
	
	
	for (int i =0; i < y.size(); i ++)
	{
		chiq += pow((y.at(i) - (a +b*x.at(i))), 2)/y.at(i); 
	}
	
	cout << "il chi quadro tra 320 335, 390 405 vale: " << chiq << endl; 

    double x1 = S; 
    double mu = 0; 
    double sigma = sqrt(var_F + var_Y); 

    double p_value = boost::math::cdf(boost::math::complement(dist, chiq));

    double p = p_value_gaussiano(x1, mu, sigma);

    cout << "il p-value della regressione vale: " << p_value << endl; 
 
    double tempo = 0; 

    cout << "Inserisci tempo di presa del file (min): "; 
    cin >> tempo; 

    out << tempo  << "  " << p_value << "\n"; 

    return 0;
}