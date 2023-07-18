#include <QObject>
#include <drumstick/qsmf.h>

struct MidiEvent
{
    long absTime;
    int  status;
    int  data1;
    int  data2;
};

class Sequence : public QList<MidiEvent>
{
public:
    virtual ~Sequence();
    void sort();
};

class SMFBuilder : public QObject
{
    Q_OBJECT
public:
    SMFBuilder();
    void run(QString fileName);
    void generate();

public slots:
    void errorHandler(const QString& errorStr);
    void trackHandler(int track);

private:
    drumstick::QSmf *m_engine;
    Sequence m_sequence;
};
