/****************************************************************************
** Meta object code from reading C++ file 'mainwindow.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.5.1)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../Finance_Lab/mainwindow.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'mainwindow.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.5.1. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_MainWindow_t {
    QByteArrayData data[16];
    char stringdata0[475];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_MainWindow_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_MainWindow_t qt_meta_stringdata_MainWindow = {
    {
QT_MOC_LITERAL(0, 0, 10), // "MainWindow"
QT_MOC_LITERAL(1, 11, 28), // "european_call_button_pressed"
QT_MOC_LITERAL(2, 40, 0), // ""
QT_MOC_LITERAL(3, 41, 27), // "asian_option_button_pressed"
QT_MOC_LITERAL(4, 69, 30), // "lookback_option_button_pressed"
QT_MOC_LITERAL(5, 100, 28), // "down_out_call_button_pressed"
QT_MOC_LITERAL(6, 129, 33), // "lookback_integrand_button_pre..."
QT_MOC_LITERAL(7, 163, 38), // "down_out_call_integrand_butto..."
QT_MOC_LITERAL(8, 202, 41), // "arithmetic_asian_integrand_bu..."
QT_MOC_LITERAL(9, 244, 37), // "uniform_random_numbers_button..."
QT_MOC_LITERAL(10, 282, 30), // "halton_sequence_button_pressed"
QT_MOC_LITERAL(11, 313, 25), // "box_muller_button_pressed"
QT_MOC_LITERAL(12, 339, 31), // "trapezoidal_rule_button_pressed"
QT_MOC_LITERAL(13, 371, 29), // "gauss_legendre_button_pressed"
QT_MOC_LITERAL(14, 401, 30), // "clenshaw_curtis_button_pressed"
QT_MOC_LITERAL(15, 432, 42) // "wiener_and_asset_simulation_b..."

    },
    "MainWindow\0european_call_button_pressed\0"
    "\0asian_option_button_pressed\0"
    "lookback_option_button_pressed\0"
    "down_out_call_button_pressed\0"
    "lookback_integrand_button_pressed\0"
    "down_out_call_integrand_button_pressed\0"
    "arithmetic_asian_integrand_button_pressed\0"
    "uniform_random_numbers_button_pressed\0"
    "halton_sequence_button_pressed\0"
    "box_muller_button_pressed\0"
    "trapezoidal_rule_button_pressed\0"
    "gauss_legendre_button_pressed\0"
    "clenshaw_curtis_button_pressed\0"
    "wiener_and_asset_simulation_button_pressed"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_MainWindow[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
      14,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    0,   84,    2, 0x08 /* Private */,
       3,    0,   85,    2, 0x08 /* Private */,
       4,    0,   86,    2, 0x08 /* Private */,
       5,    0,   87,    2, 0x08 /* Private */,
       6,    0,   88,    2, 0x08 /* Private */,
       7,    0,   89,    2, 0x08 /* Private */,
       8,    0,   90,    2, 0x08 /* Private */,
       9,    0,   91,    2, 0x08 /* Private */,
      10,    0,   92,    2, 0x08 /* Private */,
      11,    0,   93,    2, 0x08 /* Private */,
      12,    0,   94,    2, 0x08 /* Private */,
      13,    0,   95,    2, 0x08 /* Private */,
      14,    0,   96,    2, 0x08 /* Private */,
      15,    0,   97,    2, 0x08 /* Private */,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,

       0        // eod
};

void MainWindow::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        MainWindow *_t = static_cast<MainWindow *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->european_call_button_pressed(); break;
        case 1: _t->asian_option_button_pressed(); break;
        case 2: _t->lookback_option_button_pressed(); break;
        case 3: _t->down_out_call_button_pressed(); break;
        case 4: _t->lookback_integrand_button_pressed(); break;
        case 5: _t->down_out_call_integrand_button_pressed(); break;
        case 6: _t->arithmetic_asian_integrand_button_pressed(); break;
        case 7: _t->uniform_random_numbers_button_pressed(); break;
        case 8: _t->halton_sequence_button_pressed(); break;
        case 9: _t->box_muller_button_pressed(); break;
        case 10: _t->trapezoidal_rule_button_pressed(); break;
        case 11: _t->gauss_legendre_button_pressed(); break;
        case 12: _t->clenshaw_curtis_button_pressed(); break;
        case 13: _t->wiener_and_asset_simulation_button_pressed(); break;
        default: ;
        }
    }
    Q_UNUSED(_a);
}

const QMetaObject MainWindow::staticMetaObject = {
    { &QMainWindow::staticMetaObject, qt_meta_stringdata_MainWindow.data,
      qt_meta_data_MainWindow,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR}
};


const QMetaObject *MainWindow::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *MainWindow::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_MainWindow.stringdata0))
        return static_cast<void*>(const_cast< MainWindow*>(this));
    return QMainWindow::qt_metacast(_clname);
}

int MainWindow::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QMainWindow::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 14)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 14;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 14)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 14;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
