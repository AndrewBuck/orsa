#include <orsaQt/debug.h>

#include <orsaQt/event.h>

#include <gmpxx.h>

#include <QCoreApplication>
#include <QPushButton>
#include <QVBoxLayout>

orsaQt::DebugWidget::DebugWidget() : QWidget() {
  
  setWindowTitle("console");
  
  QVBoxLayout * vlay = new QVBoxLayout(this); 
  //
  vlay->setMargin(3);
  
  _te = new QTextEdit(this);
  _te->setTextInteractionFlags(Qt::TextSelectableByMouse | Qt::TextSelectableByKeyboard);
  //
  {
    // fixed font
    QFont font = _te->currentFont();
    font.setFamily("fixed");
    font.setFixedPitch(true);
    _te->setCurrentFont(font);
  }
  //
  vlay->addWidget(_te);
  
  QHBoxLayout * hbl = new QHBoxLayout;
  //
  hbl->addStretch();
  //
  {
    QPushButton * clearPB = new QPushButton(tr("clear"),this);
    connect(clearPB,SIGNAL(clicked()),_te,SLOT(clear()));
    hbl->addWidget(clearPB);
  }
  //
  /* {
     QPushButton * closePB = new QPushButton("close",this);
     connect(closePB,SIGNAL(clicked()),this,SLOT(hide()));
     hbl->addWidget(closePB);
     }
  */
  //
  vlay->addLayout(hbl);
  
  resize(800,minimumSizeHint().height());
}

void orsaQt::DebugWidget::append(const QString & s) {
  show();
  _te->append(s);
}

void orsaQt::DebugWidget::customEvent(QEvent * e) {
  if (e->type() == __ORSAQT_EVENT_DEBUG__) {
    DebugEvent * de = (DebugEvent*)(e);
    append(de->debugMessage());
  }
}

//

orsaQt::Debug::Debug() : orsa::Debug() {
  _dw = new DebugWidget;
}

orsaQt::Debug::~Debug() {
  delete _dw;
  _dw = 0;
  // _instance = 0;
}

orsa::Debug * orsaQt::Debug::instance() {
  
  if (_instance == 0) {
    Debug::_instance = new orsaQt::Debug;
  }
  
  return Debug::_instance;
}

void orsaQt::Debug::printOut(const std::string & s) {
  
  QCoreApplication::postEvent(_dw, 
			      new orsaQt::DebugEvent(s.c_str()));
  
}
