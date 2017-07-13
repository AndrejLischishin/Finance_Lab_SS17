#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "math.h"

#include "header_files/random_functions.hpp"
#include "header_files/simulation_functions.hpp"
#include "header_files/exotic_options.hpp"
#include "header_files/multivariate_integration.hpp"

#include <string>

#include <QDebug>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    connect(ui->european_call_button, SIGNAL(released()), this, SLOT(european_call_button_pressed()));

    connect(ui->asian_button, SIGNAL(released()), this, SLOT(asian_option_button_pressed()));

    connect(ui->lookback_button, SIGNAL(released()), this, SLOT(lookback_option_button_pressed()));

    connect(ui->down_out_call_button, SIGNAL(released()), this, SLOT(down_out_call_button_pressed()));

    connect(ui->lookback_integrand_button, SIGNAL(released()), this, SLOT(lookback_integrand_button_pressed()));

    connect(ui->down_out_call_integrand_button, SIGNAL(released()), this, SLOT(down_out_call_integrand_button_pressed()));

    connect(ui->arithmetic_asian_integrand_button, SIGNAL(released()), this, SLOT(arithmetic_asian_integrand_button_pressed()));

    connect(ui->uniform_random_numbers_button, SIGNAL(released()), this, SLOT(uniform_random_numbers_button_pressed()));

    connect(ui->halton_sequence_button, SIGNAL(released()), this, SLOT(halton_sequence_button_pressed()));

    connect(ui->box_muller_button, SIGNAL(released()), this, SLOT(box_muller_button_pressed()));

    connect(ui->trapezoidal_rule_button, SIGNAL(released()), this, SLOT(trapezoidal_rule_button_pressed()));

    connect(ui->gauss_legendre_button, SIGNAL(released()), this, SLOT(gauss_legendre_button_pressed()));

    connect(ui->clenshaw_curtis_button, SIGNAL(released()), this, SLOT(clenshaw_curtis_button_pressed()));

    connect(ui->wiener_and_asset_simulation_button, SIGNAL(released()), this, SLOT(wiener_and_asset_simulation_button_pressed()));
}

MainWindow::~MainWindow()
{
    delete ui;
}

double call_option_exact_expected_value(double s0, double mu, double T, double sigma, double K)
{
    double chi = (log(K/s0)-((mu-sigma*sigma*0.5)*T))/(sigma*sqrt(T));
    return s0*exp(mu*T)*normal_cdf(sigma*sqrt(T)-chi)-K*normal_cdf(-chi);
}

void MainWindow::european_call_button_pressed()
{
    double s0 = (ui->input_s0->text()).toDouble();
    double mu = (ui->input_mu->text()).toDouble();
    double sigma = (ui->input_sigma->text()).toDouble();
    double T = (ui->input_T->text()).toDouble();
    double K = (ui->input_K->text()).toDouble();

    double exact_result = call_option_exact_expected_value(s0, mu, T, sigma, K);

    QString exact_result_label = QString::number(exact_result);

    ui->fair_price->setText("Fair price for European Call Option: "+exact_result_label);
}

void MainWindow::asian_option_button_pressed()
{
    double s0 = (ui->input_s0->text()).toDouble();
    double r = (ui->input_mu->text()).toDouble();
    double sigma = (ui->input_sigma->text()).toDouble();
    double T = (ui->input_T->text()).toDouble();
    double K = (ui->input_K->text()).toDouble();

    double exact_result = continuous_geometric_average_exact(s0, r, T, K, sigma);

    QString exact_result_label = QString::number(exact_result);

    ui->fair_price->setText("Fair price for Continuous Geometric Asian Option: "+exact_result_label);
}

void MainWindow::lookback_option_button_pressed()
{

}

void MainWindow::down_out_call_button_pressed()
{
    double s0 = (ui->input_s0->text()).toDouble();
    double r = (ui->input_mu->text()).toDouble();
    double sigma = (ui->input_sigma->text()).toDouble();
    double T = (ui->input_T->text()).toDouble();
    double K = (ui->input_K->text()).toDouble();
    double B = (ui->input_B->text()).toDouble();

    double exact_result = black_scholes_down_out_call(s0, K, T, sigma, r, B);

    QString exact_result_label = QString::number(exact_result);

    ui->fair_price->setText("Fair price for Down-Out Call Option: "+exact_result_label);
}

void MainWindow::lookback_integrand_button_pressed()
{
    double s0 = (ui->input_s0->text()).toDouble();
    double r = (ui->input_mu->text()).toDouble();
    double sigma = (ui->input_sigma->text()).toDouble();
    double T = (ui->input_T->text()).toDouble();
    double K = (ui->input_K->text()).toDouble();
    double payoff = 0.0;

    std::vector<double>* x;
    x = new std::vector<double>(2);

    std::ofstream myfile;

    myfile.open("output/plot_integrand_lookback.txt",std::ios::trunc);
    if (!myfile.is_open()) {
        std::cout<<"Error opening the file"<<std::endl;
    }

    for(int i=1; i<100; i++)
    {
        for(int j=1; j<100; j++)
        {
            x->clear();
            x->push_back(i*T/100.);
            x->push_back(j*T/100.);
            payoff = payoff_discrete_lookback(x, s0, r, T, 2, K, sigma);
            myfile << (*x)[0] << "	" << (*x)[1] << "	" << payoff << std::endl;
        }
    }

    free(x);

    myfile.close();

    FILE *stream;

    stream=popen("gnuplot", "w");
    fprintf(stream, "set terminal png size 950,500\n");
    fprintf(stream, "set output 'images/lookback_integrand.png'\n");
    fprintf(stream, "set key off\n");
    fprintf(stream, "set grid\n");
    fprintf(stream, "set title \"Payoff of discrete Lookback option\"\n");
    fprintf(stream, "splot 'output/plot_integrand_lookback.txt' using 1:2:3 with points palette pointsize 0.5 pointtype 7\n");
    fclose(stream);

    QPixmap pix("images/lookback_integrand.png");
    ui->pic_label->setPixmap(pix);
}

void MainWindow::down_out_call_integrand_button_pressed()
{
    double s0 = (ui->input_s0->text()).toDouble();
    double r = (ui->input_mu->text()).toDouble();
    double sigma = (ui->input_sigma->text()).toDouble();
    double T = (ui->input_T->text()).toDouble();
    double K = (ui->input_K->text()).toDouble();
    double B = (ui->input_B->text()).toDouble();
    double payoff = 0.0;

    std::vector<double>* x;
    x = new std::vector<double>(2);

    std::ofstream myfile;

    myfile.open("output/plot_integrand_down_out_call.txt",std::ios::trunc);
    if (!myfile.is_open()) {
        std::cout<<"Error opening the file"<<std::endl;
    }

    for(int i=1; i<100; i++)
    {
        for(int j=1; j<100; j++)
        {
            x->clear();
            x->push_back(i*T/100.);
            x->push_back(j*T/100.);
            payoff = payoff_discrete_down_out_call(x, s0, r, T, 2, K, sigma, B);
            myfile << (*x)[0] << "	" << (*x)[1] << "	" << payoff << std::endl;
        }
    }

    free(x);

    myfile.close();

    FILE *stream;

    stream=popen("gnuplot", "w");
    fprintf(stream, "set terminal png size 950,500\n");
    fprintf(stream, "set output 'images/down_out_call_integrand.png'\n");
    fprintf(stream, "set key off\n");
    fprintf(stream, "set grid\n");
    fprintf(stream, "set title \"Payoff of discrete Down-Out Call option\"\n");
    fprintf(stream, "splot 'output/plot_integrand_down_out_call.txt' using 1:2:3 with points palette pointsize 0.5 pointtype 7\n");
    fclose(stream);

    QPixmap pix("images/down_out_call_integrand.png");
    ui->pic_label->setPixmap(pix);
}

void MainWindow::arithmetic_asian_integrand_button_pressed()
{
    double s0 = (ui->input_s0->text()).toDouble();
    double r = (ui->input_mu->text()).toDouble();
    double sigma = (ui->input_sigma->text()).toDouble();
    double T = (ui->input_T->text()).toDouble();
    double K = (ui->input_K->text()).toDouble();
    double payoff = 0.0;

    std::vector<double> x;

    std::ofstream myfile;

    myfile.open("output/plot_payoff_discrete_arithmetic_average.txt",std::ios::trunc);
    if (!myfile.is_open()) {
        std::cout<<"Error opening the file"<<std::endl;
    }

    for(int i=1; i<100; i++)
    {
        for(int j=1; j<100; j++)
        {
            x.clear();
            x.push_back(i*T/100.);
            x.push_back(j*T/100.);
            payoff = payoff_discrete_arithmetic_average(x, s0, r, T, 2, K, sigma);
            myfile << x[0] << "	" << x[1] << "	" << payoff << std::endl;
        }
    }

    myfile.close();

    FILE *stream;

    stream=popen("gnuplot", "w");
    fprintf(stream, "set terminal png size 950,500\n");
    fprintf(stream, "set output 'images/arithmetic_asian_integrand.png'\n");
    fprintf(stream, "set key off\n");
    fprintf(stream, "set grid\n");
    fprintf(stream, "set title \"Payoff of discrete arithmetic Asian option\"\n");
    fprintf(stream, "splot 'output/plot_payoff_discrete_arithmetic_average.txt' using 1:2:3 with points palette pointsize 0.5 pointtype 7\n");
    fclose(stream);

    QPixmap pix("images/arithmetic_asian_integrand.png");
    ui->pic_label->setPixmap(pix);
}

void MainWindow::uniform_random_numbers_button_pressed()
{
    std::ofstream myfile;

    myfile.open("output/uniform_random_numbers.txt",std::ios::trunc);
    if (!myfile.is_open()) {
        std::cout<<"Error opening the file"<<std::endl;
    }

    int num = (ui->input_number_of_random_points->text()).toInt();

    for(int i=0; i<num; i++)
    {
        for(int j=0; j<2; j++)
        {
            myfile << random_number_01() << "	";
        }
        myfile << std::endl;
    }

    myfile.close();

    FILE *stream;

    stream=popen("gnuplot", "w");
    fprintf(stream, "set terminal png size 950,500\n");
    fprintf(stream, "set output 'images/uniform_random_numbers.png'\n");
    fprintf(stream, "set title \"Uniform random numbers\"\n");
    fprintf(stream, "plot 'output/uniform_random_numbers.txt' notitle\n");
    fclose(stream);

    QPixmap pix("images/uniform_random_numbers.png");
    ui->random_functions_pic_label->setPixmap(pix);
}

void MainWindow::halton_sequence_button_pressed()
{
    std::ofstream myfile;

    myfile.open("output/halton_sequence.txt",std::ios::trunc);
    if (!myfile.is_open()) {
        std::cout<<"Error opening the file"<<std::endl;
    }

    int num = (ui->input_number_of_random_points->text()).toInt();

    std::vector<std::vector<double>> halton_sequence;
    halton_sequence = d_dimensional_halton_sequence(2, num);
    for(int i=0; i<num; i++)
    {
        for(int j=0; j<2; j++)
        {
            myfile << halton_sequence[i][j] << "	";
        }
        myfile << std::endl;
    }

    myfile.close();

    FILE *stream;

    stream=popen("gnuplot", "w");
    fprintf(stream, "set terminal png size 950,500\n");
    fprintf(stream, "set output 'images/halton_sequence.png'\n");
    fprintf(stream, "set title \"Halton sequence\"\n");
    fprintf(stream, "plot 'output/halton_sequence.txt' notitle\n");
    fclose(stream);

    QPixmap pix("images/halton_sequence.png");
    ui->random_functions_pic_label->setPixmap(pix);
}

void MainWindow::box_muller_button_pressed()
{
    gsl_rng* r;

    //seeding
    unsigned long seed = time(NULL);
    //memory allocation
    r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r, seed);

    std::ofstream myfile;

    int num_sampl = (ui->input_number_of_random_points->text()).toInt();

    //opens a file and checks if it was successfully
    //when file will be opened previous content will be deleted
    myfile.open("output/box_muller.txt", std::ios::trunc);
    if (!myfile.is_open()) {
        std::cout << "Error opening the file" << std::endl;
    }

    // output of the mueller_box_algo in order to write it to the file
    std::vector<double>* output;
    for (int i = 0; i<num_sampl; i++) {
        output = box_muller_algo(r);
        myfile << (*output)[0] << " " << (*output)[1] << std::endl;
    }

    //closes file
    myfile.close();

    FILE *stream;

    stream=popen("gnuplot", "w");
    fprintf(stream, "set terminal png size 950,500\n");
    fprintf(stream, "set output 'images/box_muller.png'\n");
    fprintf(stream, "set title \"Box muller\"\n");
    fprintf(stream, "plot 'output/box_muller.txt' notitle\n");
    fclose(stream);

    QPixmap pix("images/box_muller.png");
    ui->random_functions_pic_label->setPixmap(pix);
}

void MainWindow::trapezoidal_rule_button_pressed()
{
    int d = 2;
    int l = (ui->level_input->text()).toInt();

    std::ofstream myfile;

    std::vector<int> Nl;
    std::vector<double>* nodes;
    std::vector<double>* weights;
    std::vector<std::vector<double>> nodes_temp;
    std::vector<int> ids;
    nodes = new std::vector<double>;
    weights = new std::vector<double>;
    trap_rule(nodes, weights, l);

    for(int i=0; i<d; i++)
    {
        Nl.push_back((int)pow(2,l)-1);
        std::vector<double> row;
        nodes_temp.push_back(row);
        ids.push_back(0);
    }

    for(int i=0; i<d; i++)
    {
        for(int j=0; j<Nl[i]; j++)
        {
            nodes_temp[i].push_back((*nodes)[j]);
        }
    }

    myfile.open("output/quadrature_points_trapezoidal_rule.txt",std::ios::trunc);
    if (!myfile.is_open()) {
        std::cout<<"Error opening the file"<<std::endl;
    }
    write_quadrature_points_to_file(myfile, 0, nodes_temp, d, Nl, ids);
    myfile.close();

    FILE *stream;

    stream=popen("gnuplot", "w");
    fprintf(stream, "set terminal png size 950,500\n");
    fprintf(stream, "set output 'images/quadrature_points_trapezoidal_rule.png'\n");
    fprintf(stream, "set title \"Trapezoidal rule quadrature points\"\n");
    fprintf(stream, "plot 'output/quadrature_points_trapezoidal_rule.txt' notitle\n");
    fclose(stream);

    QPixmap pix("images/quadrature_points_trapezoidal_rule.png");
    ui->quadrature_points_pic_label->setPixmap(pix);
}

void MainWindow::gauss_legendre_button_pressed()
{
    int d = 2;
    int l = (ui->level_input->text()).toInt();

    std::ofstream myfile;

    std::vector<int> Nl;
    std::vector<double>* nodes;
    std::vector<double>* weights;
    std::vector<std::vector<double>> nodes_temp;
    std::vector<int> ids;
    nodes = new std::vector<double>;
    weights = new std::vector<double>;
    gauss_legendre(nodes, weights, l);

    for(int i=0; i<d; i++)
    {
        Nl.push_back((int)pow(2,l)-1);
        std::vector<double> row;
        nodes_temp.push_back(row);
        ids.push_back(0);
    }

    for(int i=0; i<d; i++)
    {
        for(int j=0; j<Nl[i]; j++)
        {
            nodes_temp[i].push_back((*nodes)[j]);
        }
    }

    myfile.open("output/quadrature_points_gauss_legendre.txt",std::ios::trunc);
    if (!myfile.is_open()) {
        std::cout<<"Error opening the file"<<std::endl;
    }
    write_quadrature_points_to_file(myfile, 0, nodes_temp, d, Nl, ids);
    myfile.close();

    FILE *stream;

    stream=popen("gnuplot", "w");
    fprintf(stream, "set terminal png size 950,500\n");
    fprintf(stream, "set output 'images/quadrature_points_gauss_legendre.png'\n");
    fprintf(stream, "set title \"Gauss Legendre quadrature points\"\n");
    fprintf(stream, "plot 'output/quadrature_points_gauss_legendre.txt' notitle\n");
    fclose(stream);

    QPixmap pix("images/quadrature_points_gauss_legendre.png");
    ui->quadrature_points_pic_label->setPixmap(pix);
}

void MainWindow::clenshaw_curtis_button_pressed()
{
    int d = 2;
    int l = (ui->level_input->text()).toInt();

    std::ofstream myfile;

    std::vector<int> Nl;
    std::vector<double>* nodes;
    std::vector<double>* weights;
    std::vector<std::vector<double>> nodes_temp;
    std::vector<int> ids;
    nodes = new std::vector<double>;
    weights = new std::vector<double>;
    clenshaw_curtis(nodes, weights, l);

    for(int i=0; i<d; i++)
    {
        Nl.push_back((int)pow(2,l)-1);
        std::vector<double> row;
        nodes_temp.push_back(row);
        ids.push_back(0);
    }

    for(int i=0; i<d; i++)
    {
        for(int j=0; j<Nl[i]; j++)
        {
            nodes_temp[i].push_back((*nodes)[j]);
        }
    }

    myfile.open("output/quadrature_points_clenshaw_curtis.txt",std::ios::trunc);
    if (!myfile.is_open()) {
        std::cout<<"Error opening the file"<<std::endl;
    }
    write_quadrature_points_to_file(myfile, 0, nodes_temp, d, Nl, ids);
    myfile.close();

    FILE *stream;

    stream=popen("gnuplot", "w");
    fprintf(stream, "set terminal png size 950,500\n");
    fprintf(stream, "set output 'images/quadrature_points_clenshaw_curtis.png'\n");
    fprintf(stream, "set title \"Clenshaw Curtis quadrature points\"\n");
    fprintf(stream, "plot 'output/quadrature_points_clenshaw_curtis.txt' notitle\n");
    fclose(stream);

    QPixmap pix("images/quadrature_points_clenshaw_curtis.png");
    ui->quadrature_points_pic_label->setPixmap(pix);
}

void MainWindow::wiener_and_asset_simulation_button_pressed()
{
    gsl_rng* r;

    //seeding
    unsigned long seed = time(NULL);
    //memory allocation
    r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r, seed);

    int number_of_simulations = (ui->input_number_of_simulations->text()).toInt();
    double delta_t = (ui->input_delta_t->text()).toDouble();

    double s0 = (ui->input_s0->text()).toDouble();
    double mu = (ui->input_mu->text()).toDouble();
    double sigma = (ui->input_sigma->text()).toDouble();
    double T = (ui->input_T->text()).toDouble();

    std::ofstream wiener_file;
    std::ofstream asset_file;

    wiener_file.open("output/wiener_process.txt",std::ios::trunc);
    if (!wiener_file.is_open()) {
        std::cout<<"Error opening the file"<<std::endl;
    }
    int M = (int)(T/delta_t);

    std::vector<std::vector<double>*> wiener_all;

    for(int j=0; j<number_of_simulations; j++)
    {
        std::vector<double>* w = wiener_process(r, T, delta_t);
        wiener_all.push_back(w);
    }

    for(int i=0; i<=M; i++)
    {
        wiener_file<<i*delta_t<<"   ";
        for(int k=0; k<number_of_simulations; k++)
        {
            wiener_file<<(*wiener_all[k])[i]<<" ";
        }
        wiener_file<<std::endl;
    }

    wiener_file.close();

    FILE *stream;

    stream=popen("gnuplot", "w");
    fprintf(stream, "set terminal png size 950,250\n");
    fprintf(stream, "set output 'images/wiener_process.png'\n");
    fprintf(stream, "set title \"Wiener process\"\n");
    std::string temp = "plot ";
    for(int i=0; i<number_of_simulations; i++)
    {
        temp += "\"output/wiener_process.txt\" using 1:";
        temp += std::to_string(i+2);
        temp += " with linespoints ls 1 notitle";
        if(i!=number_of_simulations-1)
            temp += ", ";
    }
    fprintf(stream, (temp).c_str());
    fclose(stream);

    QPixmap pix("images/wiener_process.png");
    ui->wiener_process_pic_label->setPixmap(pix);

    asset_file.open("output/correspond_asset_prices.txt",std::ios::trunc);
    if (!asset_file.is_open()) {
        std::cout<<"Error opening the file"<<std::endl;
    }

    std::vector<std::vector<double>*> asset_prices_all;

    for(int j=0; j<number_of_simulations; j++)
    {
        std::vector<double>* s = brownian_motion(r, T, delta_t, wiener_all[j], s0, mu, sigma);
        asset_prices_all.push_back(s);
    }

    for(int i=0; i<=M; i++)
    {
        asset_file<<i*delta_t<<" ";
        for(int k=0; k<number_of_simulations; k++)
        {
            asset_file<<(*asset_prices_all[k])[i]<<" ";
        }
        asset_file<<std::endl;
    }

    asset_file.close();

    stream=popen("gnuplot", "w");
    fprintf(stream, "set terminal png size 950,250\n");
    fprintf(stream, "set output 'images/correspond_asset_prices.png'\n");
    fprintf(stream, "set title \"Asset prices\"\n");
    temp = "plot ";
    for(int i=0; i<number_of_simulations; i++)
    {
        temp += "\"output/correspond_asset_prices.txt\" using 1:";
        temp += std::to_string(i+2);
        temp += " with linespoints ls 1 notitle";
        if(i!=number_of_simulations-1)
            temp += ", ";
    }
    fprintf(stream, (temp).c_str());
    fclose(stream);

    QPixmap pix2("images/correspond_asset_prices.png");
    ui->asset_price_pic_label->setPixmap(pix2);
}
