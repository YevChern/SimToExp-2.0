/****************************************************************************
** Meta object code from reading C++ file 'ExpManager.h'
**
** Created by: The Qt Meta Object Compiler version 63 (Qt 4.8.7)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../ExpManager.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'ExpManager.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 63
#error "This file was generated using the moc from 4.8.7. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_ExpManager[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
       3,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       2,       // signalCount

 // signals: signature, parameters, type, tag, flags
      12,   11,   11,   11, 0x05,
      29,   11,   11,   11, 0x05,

 // slots: signature, parameters, type, tag, flags
      49,   11,   11,   11, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_ExpManager[] = {
    "ExpManager\0\0xrayDataScaled()\0"
    "neutronDataScaled()\0resetAllData()\0"
};

void ExpManager::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Q_ASSERT(staticMetaObject.cast(_o));
        ExpManager *_t = static_cast<ExpManager *>(_o);
        switch (_id) {
        case 0: _t->xrayDataScaled(); break;
        case 1: _t->neutronDataScaled(); break;
        case 2: _t->resetAllData(); break;
        default: ;
        }
    }
    Q_UNUSED(_a);
}

const QMetaObjectExtraData ExpManager::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject ExpManager::staticMetaObject = {
    { &QWidget::staticMetaObject, qt_meta_stringdata_ExpManager,
      qt_meta_data_ExpManager, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &ExpManager::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *ExpManager::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *ExpManager::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_ExpManager))
        return static_cast<void*>(const_cast< ExpManager*>(this));
    return QWidget::qt_metacast(_clname);
}

int ExpManager::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 3)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 3;
    }
    return _id;
}

// SIGNAL 0
void ExpManager::xrayDataScaled()
{
    QMetaObject::activate(this, &staticMetaObject, 0, 0);
}

// SIGNAL 1
void ExpManager::neutronDataScaled()
{
    QMetaObject::activate(this, &staticMetaObject, 1, 0);
}
QT_END_MOC_NAMESPACE
