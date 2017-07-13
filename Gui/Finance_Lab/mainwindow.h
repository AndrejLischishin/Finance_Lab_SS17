#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private:
    Ui::MainWindow *ui;

private slots:
    void european_call_button_pressed();

private slots:
    void asian_option_button_pressed();

private slots:
    void lookback_option_button_pressed();

private slots:
    void down_out_call_button_pressed();

private slots:
    void lookback_integrand_button_pressed();

private slots:
    void down_out_call_integrand_button_pressed();

private slots:
    void arithmetic_asian_integrand_button_pressed();

private slots:
    void uniform_random_numbers_button_pressed();

private slots:
    void halton_sequence_button_pressed();

private slots:
    void box_muller_button_pressed();

private slots:
    void trapezoidal_rule_button_pressed();

private slots:
    void gauss_legendre_button_pressed();

private slots:
    void clenshaw_curtis_button_pressed();

private slots:
    void wiener_and_asset_simulation_button_pressed();
};

#endif // MAINWINDOW_H
