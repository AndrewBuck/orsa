#ifndef _ORSAQT_DEBUG_
#define _ORSAQT_DEBUG_

#include <QTextEdit>
#include <QWidget>

#include <orsa/debug.h>

namespace orsaQt {
  
  class DebugWidget : public QWidget {
    
    Q_OBJECT;
    
  public:	
    DebugWidget();
    
  public:
    virtual ~DebugWidget() {
      ORSA_DEBUG("~DebugWidget() called...");
    }
    
  public slots:
    void append(const QString &);
    
  protected:
    void customEvent(QEvent *);
    
  private:
    QTextEdit * _te;
  };
  
  class Debug : public orsa::Debug {
    
  protected:
    Debug();
    
  public:
    virtual ~Debug();
    
  public:
    static orsa::Debug * instance();
    
  protected:
    void printOut(const std::string &);
    
  private:
    DebugWidget * _dw;
  };
  
} // namespace orsaQt

#endif // _ORSAQT_DEBUG_
