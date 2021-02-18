/****************************************************************************
** Meta object code from reading C++ file 'TrajectoryInfoWindow.h'
**
** Created by: The Qt Meta Object Compiler version 63 (Qt 4.8.7)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../TrajectoryInfoWindow.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'TrajectoryInfoWindow.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 63
#error "This file was generated using the moc from 4.8.7. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_TrajectoryInfoWindow[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
      12,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       2,       // signalCount

 // signals: signature, parameters, type, tag, flags
      22,   21,   21,   21, 0x05,
      57,   21,   21,   21, 0x05,

 // slots: signature, parameters, type, tag, flags
      89,   21,   21,   21, 0x08,
     113,   21,   21,   21, 0x08,
     140,   21,   21,   21, 0x08,
     158,   21,   21,   21, 0x08,
     178,   21,   21,   21, 0x08,
     207,   21,   21,   21, 0x08,
     236,   21,   21,   21, 0x08,
     253,   21,   21,   21, 0x08,
     270,   21,   21,   21, 0x08,
     290,   21,   21,   21, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_TrajectoryInfoWindow[] = {
    "TrajectoryInfoWindow\0\0"
    "dataEntered(TrajectoryInfoWindow*)\0"
    "canceled(TrajectoryInfoWindow*)\0"
    "updateFilePath(QString)\0"
    "updateOptFilePath(QString)\0openFileBrowser()\0"
    "updateTrajInfo(int)\0normalLineEditWritable(bool)\0"
    "planarLineEditWritable(bool)\0"
    "calcNormal(bool)\0calcPlanar(bool)\0"
    "dataEntryComplete()\0cancel()\0"
};

void TrajectoryInfoWindow::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Q_ASSERT(staticMetaObject.cast(_o));
        TrajectoryInfoWindow *_t = static_cast<TrajectoryInfoWindow *>(_o);
        switch (_id) {
        case 0: _t->dataEntered((*reinterpret_cast< TrajectoryInfoWindow*(*)>(_a[1]))); break;
        case 1: _t->canceled((*reinterpret_cast< TrajectoryInfoWindow*(*)>(_a[1]))); break;
        case 2: _t->updateFilePath((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 3: _t->updateOptFilePath((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 4: _t->openFileBrowser(); break;
        case 5: _t->updateTrajInfo((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 6: _t->normalLineEditWritable((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 7: _t->planarLineEditWritable((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 8: _t->calcNormal((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 9: _t->calcPlanar((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 10: _t->dataEntryComplete(); break;
        case 11: _t->cancel(); break;
        default: ;
        }
    }
}

const QMetaObjectExtraData TrajectoryInfoWindow::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject TrajectoryInfoWindow::staticMetaObject = {
    { &QWidget::staticMetaObject, qt_meta_stringdata_TrajectoryInfoWindow,
      qt_meta_data_TrajectoryInfoWindow, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &TrajectoryInfoWindow::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *TrajectoryInfoWindow::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *TrajectoryInfoWindow::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_TrajectoryInfoWindow))
        return static_cast<void*>(const_cast< TrajectoryInfoWindow*>(this));
    return QWidget::qt_metacast(_clname);
}

int TrajectoryInfoWindow::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 12)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 12;
    }
    return _id;
}

// SIGNAL 0
void TrajectoryInfoWindow::dataEntered(TrajectoryInfoWindow * _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}

// SIGNAL 1
void TrajectoryInfoWindow::canceled(TrajectoryInfoWindow * _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 1, _a);
}
QT_END_MOC_NAMESPACE
