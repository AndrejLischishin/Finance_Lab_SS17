/********************************************************************************
** Form generated from reading UI file 'mainwindow.ui'
**
** Created by: Qt User Interface Compiler version 5.5.1
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_MAINWINDOW_H
#define UI_MAINWINDOW_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QSpacerItem>
#include <QtWidgets/QTabWidget>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_MainWindow
{
public:
    QWidget *centralWidget;
    QHBoxLayout *horizontalLayout;
    QTabWidget *tabWidget;
    QWidget *tab;
    QHBoxLayout *horizontalLayout_2;
    QGridLayout *gridLayout_2;
    QLineEdit *input_K;
    QLabel *label;
    QLabel *label_6;
    QLabel *label_2;
    QLabel *label_4;
    QLineEdit *input_sigma;
    QLineEdit *input_T;
    QLabel *label_3;
    QLineEdit *input_s0;
    QLineEdit *input_mu;
    QLineEdit *input_B;
    QLabel *label_5;
    QVBoxLayout *verticalLayout;
    QPushButton *european_call_button;
    QPushButton *asian_button;
    QPushButton *down_out_call_button;
    QPushButton *lookback_button;
    QLabel *fair_price;
    QWidget *tab_5;
    QVBoxLayout *verticalLayout_7;
    QHBoxLayout *horizontalLayout_9;
    QLabel *label_10;
    QLineEdit *input_delta_t;
    QLabel *label_11;
    QLineEdit *input_number_of_simulations;
    QPushButton *wiener_and_asset_simulation_button;
    QSpacerItem *horizontalSpacer_3;
    QLabel *wiener_process_pic_label;
    QLabel *asset_price_pic_label;
    QWidget *tab_2;
    QHBoxLayout *horizontalLayout_3;
    QVBoxLayout *verticalLayout_2;
    QHBoxLayout *horizontalLayout_4;
    QPushButton *arithmetic_asian_integrand_button;
    QPushButton *lookback_integrand_button;
    QPushButton *down_out_call_integrand_button;
    QLabel *pic_label;
    QWidget *tab_3;
    QVBoxLayout *verticalLayout_3;
    QVBoxLayout *verticalLayout_6;
    QHBoxLayout *horizontalLayout_5;
    QLabel *label_8;
    QLineEdit *input_number_of_random_points;
    QSpacerItem *horizontalSpacer_2;
    QHBoxLayout *horizontalLayout_8;
    QPushButton *uniform_random_numbers_button;
    QPushButton *halton_sequence_button;
    QPushButton *box_muller_button;
    QLabel *random_functions_pic_label;
    QWidget *tab_4;
    QVBoxLayout *verticalLayout_4;
    QVBoxLayout *verticalLayout_5;
    QHBoxLayout *horizontalLayout_7;
    QLabel *label_7;
    QLineEdit *level_input;
    QSpacerItem *horizontalSpacer;
    QHBoxLayout *horizontalLayout_6;
    QPushButton *trapezoidal_rule_button;
    QPushButton *gauss_legendre_button;
    QPushButton *clenshaw_curtis_button;
    QLabel *quadrature_points_pic_label;

    void setupUi(QMainWindow *MainWindow)
    {
        if (MainWindow->objectName().isEmpty())
            MainWindow->setObjectName(QStringLiteral("MainWindow"));
        MainWindow->resize(1000, 620);
        centralWidget = new QWidget(MainWindow);
        centralWidget->setObjectName(QStringLiteral("centralWidget"));
        QSizePolicy sizePolicy(QSizePolicy::Minimum, QSizePolicy::Preferred);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(centralWidget->sizePolicy().hasHeightForWidth());
        centralWidget->setSizePolicy(sizePolicy);
        centralWidget->setMinimumSize(QSize(625, 0));
        horizontalLayout = new QHBoxLayout(centralWidget);
        horizontalLayout->setSpacing(6);
        horizontalLayout->setContentsMargins(11, 11, 11, 11);
        horizontalLayout->setObjectName(QStringLiteral("horizontalLayout"));
        tabWidget = new QTabWidget(centralWidget);
        tabWidget->setObjectName(QStringLiteral("tabWidget"));
        tab = new QWidget();
        tab->setObjectName(QStringLiteral("tab"));
        horizontalLayout_2 = new QHBoxLayout(tab);
        horizontalLayout_2->setSpacing(6);
        horizontalLayout_2->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_2->setObjectName(QStringLiteral("horizontalLayout_2"));
        gridLayout_2 = new QGridLayout();
        gridLayout_2->setSpacing(6);
        gridLayout_2->setObjectName(QStringLiteral("gridLayout_2"));
        input_K = new QLineEdit(tab);
        input_K->setObjectName(QStringLiteral("input_K"));
        QSizePolicy sizePolicy1(QSizePolicy::Minimum, QSizePolicy::Fixed);
        sizePolicy1.setHorizontalStretch(0);
        sizePolicy1.setVerticalStretch(0);
        sizePolicy1.setHeightForWidth(input_K->sizePolicy().hasHeightForWidth());
        input_K->setSizePolicy(sizePolicy1);

        gridLayout_2->addWidget(input_K, 5, 1, 1, 1);

        label = new QLabel(tab);
        label->setObjectName(QStringLiteral("label"));

        gridLayout_2->addWidget(label, 0, 0, 1, 1);

        label_6 = new QLabel(tab);
        label_6->setObjectName(QStringLiteral("label_6"));

        gridLayout_2->addWidget(label_6, 5, 0, 1, 1);

        label_2 = new QLabel(tab);
        label_2->setObjectName(QStringLiteral("label_2"));

        gridLayout_2->addWidget(label_2, 1, 0, 1, 1);

        label_4 = new QLabel(tab);
        label_4->setObjectName(QStringLiteral("label_4"));

        gridLayout_2->addWidget(label_4, 4, 0, 1, 1);

        input_sigma = new QLineEdit(tab);
        input_sigma->setObjectName(QStringLiteral("input_sigma"));
        sizePolicy1.setHeightForWidth(input_sigma->sizePolicy().hasHeightForWidth());
        input_sigma->setSizePolicy(sizePolicy1);

        gridLayout_2->addWidget(input_sigma, 3, 1, 1, 1);

        input_T = new QLineEdit(tab);
        input_T->setObjectName(QStringLiteral("input_T"));
        sizePolicy1.setHeightForWidth(input_T->sizePolicy().hasHeightForWidth());
        input_T->setSizePolicy(sizePolicy1);

        gridLayout_2->addWidget(input_T, 4, 1, 1, 1);

        label_3 = new QLabel(tab);
        label_3->setObjectName(QStringLiteral("label_3"));

        gridLayout_2->addWidget(label_3, 3, 0, 1, 1);

        input_s0 = new QLineEdit(tab);
        input_s0->setObjectName(QStringLiteral("input_s0"));
        sizePolicy1.setHeightForWidth(input_s0->sizePolicy().hasHeightForWidth());
        input_s0->setSizePolicy(sizePolicy1);

        gridLayout_2->addWidget(input_s0, 0, 1, 1, 1);

        input_mu = new QLineEdit(tab);
        input_mu->setObjectName(QStringLiteral("input_mu"));
        sizePolicy1.setHeightForWidth(input_mu->sizePolicy().hasHeightForWidth());
        input_mu->setSizePolicy(sizePolicy1);

        gridLayout_2->addWidget(input_mu, 1, 1, 1, 1);

        input_B = new QLineEdit(tab);
        input_B->setObjectName(QStringLiteral("input_B"));
        sizePolicy1.setHeightForWidth(input_B->sizePolicy().hasHeightForWidth());
        input_B->setSizePolicy(sizePolicy1);

        gridLayout_2->addWidget(input_B, 6, 1, 1, 1);

        label_5 = new QLabel(tab);
        label_5->setObjectName(QStringLiteral("label_5"));

        gridLayout_2->addWidget(label_5, 6, 0, 1, 1);


        horizontalLayout_2->addLayout(gridLayout_2);

        verticalLayout = new QVBoxLayout();
        verticalLayout->setSpacing(6);
        verticalLayout->setObjectName(QStringLiteral("verticalLayout"));
        european_call_button = new QPushButton(tab);
        european_call_button->setObjectName(QStringLiteral("european_call_button"));

        verticalLayout->addWidget(european_call_button);

        asian_button = new QPushButton(tab);
        asian_button->setObjectName(QStringLiteral("asian_button"));

        verticalLayout->addWidget(asian_button);

        down_out_call_button = new QPushButton(tab);
        down_out_call_button->setObjectName(QStringLiteral("down_out_call_button"));

        verticalLayout->addWidget(down_out_call_button);

        lookback_button = new QPushButton(tab);
        lookback_button->setObjectName(QStringLiteral("lookback_button"));

        verticalLayout->addWidget(lookback_button);

        fair_price = new QLabel(tab);
        fair_price->setObjectName(QStringLiteral("fair_price"));

        verticalLayout->addWidget(fair_price);


        horizontalLayout_2->addLayout(verticalLayout);

        tabWidget->addTab(tab, QString());
        tab_5 = new QWidget();
        tab_5->setObjectName(QStringLiteral("tab_5"));
        verticalLayout_7 = new QVBoxLayout(tab_5);
        verticalLayout_7->setSpacing(6);
        verticalLayout_7->setContentsMargins(11, 11, 11, 11);
        verticalLayout_7->setObjectName(QStringLiteral("verticalLayout_7"));
        horizontalLayout_9 = new QHBoxLayout();
        horizontalLayout_9->setSpacing(6);
        horizontalLayout_9->setObjectName(QStringLiteral("horizontalLayout_9"));
        label_10 = new QLabel(tab_5);
        label_10->setObjectName(QStringLiteral("label_10"));

        horizontalLayout_9->addWidget(label_10);

        input_delta_t = new QLineEdit(tab_5);
        input_delta_t->setObjectName(QStringLiteral("input_delta_t"));
        QSizePolicy sizePolicy2(QSizePolicy::Fixed, QSizePolicy::Fixed);
        sizePolicy2.setHorizontalStretch(0);
        sizePolicy2.setVerticalStretch(0);
        sizePolicy2.setHeightForWidth(input_delta_t->sizePolicy().hasHeightForWidth());
        input_delta_t->setSizePolicy(sizePolicy2);

        horizontalLayout_9->addWidget(input_delta_t);

        label_11 = new QLabel(tab_5);
        label_11->setObjectName(QStringLiteral("label_11"));

        horizontalLayout_9->addWidget(label_11);

        input_number_of_simulations = new QLineEdit(tab_5);
        input_number_of_simulations->setObjectName(QStringLiteral("input_number_of_simulations"));
        sizePolicy2.setHeightForWidth(input_number_of_simulations->sizePolicy().hasHeightForWidth());
        input_number_of_simulations->setSizePolicy(sizePolicy2);

        horizontalLayout_9->addWidget(input_number_of_simulations);

        wiener_and_asset_simulation_button = new QPushButton(tab_5);
        wiener_and_asset_simulation_button->setObjectName(QStringLiteral("wiener_and_asset_simulation_button"));

        horizontalLayout_9->addWidget(wiener_and_asset_simulation_button);

        horizontalSpacer_3 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_9->addItem(horizontalSpacer_3);


        verticalLayout_7->addLayout(horizontalLayout_9);

        wiener_process_pic_label = new QLabel(tab_5);
        wiener_process_pic_label->setObjectName(QStringLiteral("wiener_process_pic_label"));

        verticalLayout_7->addWidget(wiener_process_pic_label);

        asset_price_pic_label = new QLabel(tab_5);
        asset_price_pic_label->setObjectName(QStringLiteral("asset_price_pic_label"));

        verticalLayout_7->addWidget(asset_price_pic_label);

        tabWidget->addTab(tab_5, QString());
        tab_2 = new QWidget();
        tab_2->setObjectName(QStringLiteral("tab_2"));
        horizontalLayout_3 = new QHBoxLayout(tab_2);
        horizontalLayout_3->setSpacing(6);
        horizontalLayout_3->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_3->setObjectName(QStringLiteral("horizontalLayout_3"));
        verticalLayout_2 = new QVBoxLayout();
        verticalLayout_2->setSpacing(6);
        verticalLayout_2->setObjectName(QStringLiteral("verticalLayout_2"));
        horizontalLayout_4 = new QHBoxLayout();
        horizontalLayout_4->setSpacing(6);
        horizontalLayout_4->setObjectName(QStringLiteral("horizontalLayout_4"));
        arithmetic_asian_integrand_button = new QPushButton(tab_2);
        arithmetic_asian_integrand_button->setObjectName(QStringLiteral("arithmetic_asian_integrand_button"));

        horizontalLayout_4->addWidget(arithmetic_asian_integrand_button);

        lookback_integrand_button = new QPushButton(tab_2);
        lookback_integrand_button->setObjectName(QStringLiteral("lookback_integrand_button"));

        horizontalLayout_4->addWidget(lookback_integrand_button);

        down_out_call_integrand_button = new QPushButton(tab_2);
        down_out_call_integrand_button->setObjectName(QStringLiteral("down_out_call_integrand_button"));

        horizontalLayout_4->addWidget(down_out_call_integrand_button);


        verticalLayout_2->addLayout(horizontalLayout_4);

        pic_label = new QLabel(tab_2);
        pic_label->setObjectName(QStringLiteral("pic_label"));

        verticalLayout_2->addWidget(pic_label);


        horizontalLayout_3->addLayout(verticalLayout_2);

        tabWidget->addTab(tab_2, QString());
        tab_3 = new QWidget();
        tab_3->setObjectName(QStringLiteral("tab_3"));
        verticalLayout_3 = new QVBoxLayout(tab_3);
        verticalLayout_3->setSpacing(6);
        verticalLayout_3->setContentsMargins(11, 11, 11, 11);
        verticalLayout_3->setObjectName(QStringLiteral("verticalLayout_3"));
        verticalLayout_6 = new QVBoxLayout();
        verticalLayout_6->setSpacing(6);
        verticalLayout_6->setObjectName(QStringLiteral("verticalLayout_6"));
        horizontalLayout_5 = new QHBoxLayout();
        horizontalLayout_5->setSpacing(6);
        horizontalLayout_5->setObjectName(QStringLiteral("horizontalLayout_5"));
        label_8 = new QLabel(tab_3);
        label_8->setObjectName(QStringLiteral("label_8"));
        QSizePolicy sizePolicy3(QSizePolicy::Fixed, QSizePolicy::Preferred);
        sizePolicy3.setHorizontalStretch(0);
        sizePolicy3.setVerticalStretch(0);
        sizePolicy3.setHeightForWidth(label_8->sizePolicy().hasHeightForWidth());
        label_8->setSizePolicy(sizePolicy3);

        horizontalLayout_5->addWidget(label_8);

        input_number_of_random_points = new QLineEdit(tab_3);
        input_number_of_random_points->setObjectName(QStringLiteral("input_number_of_random_points"));
        sizePolicy2.setHeightForWidth(input_number_of_random_points->sizePolicy().hasHeightForWidth());
        input_number_of_random_points->setSizePolicy(sizePolicy2);

        horizontalLayout_5->addWidget(input_number_of_random_points);

        horizontalSpacer_2 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_5->addItem(horizontalSpacer_2);


        verticalLayout_6->addLayout(horizontalLayout_5);

        horizontalLayout_8 = new QHBoxLayout();
        horizontalLayout_8->setSpacing(6);
        horizontalLayout_8->setObjectName(QStringLiteral("horizontalLayout_8"));
        uniform_random_numbers_button = new QPushButton(tab_3);
        uniform_random_numbers_button->setObjectName(QStringLiteral("uniform_random_numbers_button"));

        horizontalLayout_8->addWidget(uniform_random_numbers_button);

        halton_sequence_button = new QPushButton(tab_3);
        halton_sequence_button->setObjectName(QStringLiteral("halton_sequence_button"));

        horizontalLayout_8->addWidget(halton_sequence_button);

        box_muller_button = new QPushButton(tab_3);
        box_muller_button->setObjectName(QStringLiteral("box_muller_button"));

        horizontalLayout_8->addWidget(box_muller_button);


        verticalLayout_6->addLayout(horizontalLayout_8);


        verticalLayout_3->addLayout(verticalLayout_6);

        random_functions_pic_label = new QLabel(tab_3);
        random_functions_pic_label->setObjectName(QStringLiteral("random_functions_pic_label"));

        verticalLayout_3->addWidget(random_functions_pic_label);

        tabWidget->addTab(tab_3, QString());
        tab_4 = new QWidget();
        tab_4->setObjectName(QStringLiteral("tab_4"));
        verticalLayout_4 = new QVBoxLayout(tab_4);
        verticalLayout_4->setSpacing(6);
        verticalLayout_4->setContentsMargins(11, 11, 11, 11);
        verticalLayout_4->setObjectName(QStringLiteral("verticalLayout_4"));
        verticalLayout_5 = new QVBoxLayout();
        verticalLayout_5->setSpacing(6);
        verticalLayout_5->setObjectName(QStringLiteral("verticalLayout_5"));
        horizontalLayout_7 = new QHBoxLayout();
        horizontalLayout_7->setSpacing(6);
        horizontalLayout_7->setObjectName(QStringLiteral("horizontalLayout_7"));
        label_7 = new QLabel(tab_4);
        label_7->setObjectName(QStringLiteral("label_7"));
        sizePolicy3.setHeightForWidth(label_7->sizePolicy().hasHeightForWidth());
        label_7->setSizePolicy(sizePolicy3);

        horizontalLayout_7->addWidget(label_7);

        level_input = new QLineEdit(tab_4);
        level_input->setObjectName(QStringLiteral("level_input"));
        sizePolicy2.setHeightForWidth(level_input->sizePolicy().hasHeightForWidth());
        level_input->setSizePolicy(sizePolicy2);

        horizontalLayout_7->addWidget(level_input);

        horizontalSpacer = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_7->addItem(horizontalSpacer);


        verticalLayout_5->addLayout(horizontalLayout_7);

        horizontalLayout_6 = new QHBoxLayout();
        horizontalLayout_6->setSpacing(6);
        horizontalLayout_6->setObjectName(QStringLiteral("horizontalLayout_6"));
        trapezoidal_rule_button = new QPushButton(tab_4);
        trapezoidal_rule_button->setObjectName(QStringLiteral("trapezoidal_rule_button"));

        horizontalLayout_6->addWidget(trapezoidal_rule_button);

        gauss_legendre_button = new QPushButton(tab_4);
        gauss_legendre_button->setObjectName(QStringLiteral("gauss_legendre_button"));

        horizontalLayout_6->addWidget(gauss_legendre_button);

        clenshaw_curtis_button = new QPushButton(tab_4);
        clenshaw_curtis_button->setObjectName(QStringLiteral("clenshaw_curtis_button"));

        horizontalLayout_6->addWidget(clenshaw_curtis_button);


        verticalLayout_5->addLayout(horizontalLayout_6);


        verticalLayout_4->addLayout(verticalLayout_5);

        quadrature_points_pic_label = new QLabel(tab_4);
        quadrature_points_pic_label->setObjectName(QStringLiteral("quadrature_points_pic_label"));

        verticalLayout_4->addWidget(quadrature_points_pic_label);

        tabWidget->addTab(tab_4, QString());

        horizontalLayout->addWidget(tabWidget);

        MainWindow->setCentralWidget(centralWidget);

        retranslateUi(MainWindow);

        tabWidget->setCurrentIndex(1);


        QMetaObject::connectSlotsByName(MainWindow);
    } // setupUi

    void retranslateUi(QMainWindow *MainWindow)
    {
        MainWindow->setWindowTitle(QApplication::translate("MainWindow", "Finance Lab", 0));
        input_K->setText(QApplication::translate("MainWindow", "10", 0));
        label->setText(QApplication::translate("MainWindow", "S(0)", 0));
        label_6->setText(QApplication::translate("MainWindow", "K", 0));
        label_2->setText(QApplication::translate("MainWindow", "Mu/r", 0));
        label_4->setText(QApplication::translate("MainWindow", "T", 0));
        input_sigma->setText(QApplication::translate("MainWindow", "0.2", 0));
        input_T->setText(QApplication::translate("MainWindow", "1", 0));
        label_3->setText(QApplication::translate("MainWindow", "Sigma", 0));
        input_s0->setText(QApplication::translate("MainWindow", "10", 0));
        input_mu->setText(QApplication::translate("MainWindow", "0.02", 0));
        input_B->setText(QApplication::translate("MainWindow", "8", 0));
        label_5->setText(QApplication::translate("MainWindow", "B", 0));
        european_call_button->setText(QApplication::translate("MainWindow", "European Call Option", 0));
        asian_button->setText(QApplication::translate("MainWindow", "Continuous Geometric Asian Option", 0));
        down_out_call_button->setText(QApplication::translate("MainWindow", "Down-Out Call Option", 0));
        lookback_button->setText(QApplication::translate("MainWindow", "Lookback Option", 0));
        fair_price->setText(QApplication::translate("MainWindow", "Fair price:", 0));
        tabWidget->setTabText(tabWidget->indexOf(tab), QApplication::translate("MainWindow", "Parameters and Option pricing", 0));
        label_10->setText(QApplication::translate("MainWindow", "Delta t", 0));
        input_delta_t->setText(QApplication::translate("MainWindow", "0.1", 0));
        label_11->setText(QApplication::translate("MainWindow", "Number of simulations", 0));
        input_number_of_simulations->setText(QApplication::translate("MainWindow", "3", 0));
        wiener_and_asset_simulation_button->setText(QApplication::translate("MainWindow", "Simulate", 0));
        wiener_process_pic_label->setText(QString());
        asset_price_pic_label->setText(QString());
        tabWidget->setTabText(tabWidget->indexOf(tab_5), QApplication::translate("MainWindow", "Price simulations", 0));
        arithmetic_asian_integrand_button->setText(QApplication::translate("MainWindow", "Arithmetic asian option integrand", 0));
        lookback_integrand_button->setText(QApplication::translate("MainWindow", "Lookback integrand", 0));
        down_out_call_integrand_button->setText(QApplication::translate("MainWindow", "Down-Out call integrand", 0));
        pic_label->setText(QString());
        tabWidget->setTabText(tabWidget->indexOf(tab_2), QApplication::translate("MainWindow", "Payoff plots", 0));
        label_8->setText(QApplication::translate("MainWindow", "Number of random points", 0));
        input_number_of_random_points->setText(QApplication::translate("MainWindow", "100", 0));
        uniform_random_numbers_button->setText(QApplication::translate("MainWindow", "Uniform random numbers", 0));
        halton_sequence_button->setText(QApplication::translate("MainWindow", "Halton sequence", 0));
        box_muller_button->setText(QApplication::translate("MainWindow", "Box muller", 0));
        random_functions_pic_label->setText(QString());
        tabWidget->setTabText(tabWidget->indexOf(tab_3), QApplication::translate("MainWindow", "Random functions", 0));
        label_7->setText(QApplication::translate("MainWindow", "Level", 0));
        level_input->setText(QApplication::translate("MainWindow", "5", 0));
        trapezoidal_rule_button->setText(QApplication::translate("MainWindow", "Trapezoidal rule", 0));
        gauss_legendre_button->setText(QApplication::translate("MainWindow", "Gauss Legendre", 0));
        clenshaw_curtis_button->setText(QApplication::translate("MainWindow", "Clenshaw Curtis", 0));
        quadrature_points_pic_label->setText(QString());
        tabWidget->setTabText(tabWidget->indexOf(tab_4), QApplication::translate("MainWindow", "Quadrature points", 0));
    } // retranslateUi

};

namespace Ui {
    class MainWindow: public Ui_MainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_MAINWINDOW_H
