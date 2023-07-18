#include <QCoreApplication>
#include <QDebug>
#include <QTextCodec>
#include <drumstick/qsmf.h>
#include "MidiGenerator.h"

static inline bool eventLessThan(const MidiEvent& s1, const MidiEvent& s2)
{
    return s1.absTime < s2.absTime;
}

void Sequence::sort()
{
    std::stable_sort(begin(), end(), eventLessThan);
}

Sequence::~Sequence()
{
    clear();
}

SMFBuilder::SMFBuilder() : QObject()
{
    m_engine = new drumstick::QSmf(this);
    m_engine->setTextCodec(QTextCodec::codecForName("UTF-8"));
    connect(m_engine, SIGNAL(signalSMFError(const QString&)), 
            this, SLOT(errorHandler(const QString&)));
    connect(m_engine, SIGNAL(signalSMFWriteTrack(int)), 
            this, SLOT(trackHandler(int)));
}

void SMFBuilder::errorHandler(const QString& errorStr)
{
    qWarning() << errorStr;
    exit(1);
}

void SMFBuilder::trackHandler(int )
{
    // meta events
    m_engine->writeBpmTempo(0, 100); // m.m.=100 (andante)
    m_engine->writeTimeSignature(0, 4, 2, 18, 8);  // ts = 4/4
    m_engine->writeKeySignature(0, 0, major_mode); // C major
    m_engine->writeMidiEvent(0, program_chng, 0, 0);  // grand piano
    // note events
    long last_time = 0;
    for(auto ev : m_sequence)
    {
        long delta = ev.absTime - last_time;
        last_time = ev.absTime;
        m_engine->writeMidiEvent(delta,  ev.status, 0, ev.data1, ev.data2);
    }
    // final event
    m_engine->writeMetaEvent(0, end_of_track); 
}

void SMFBuilder::run(QString fileName)
{
    m_engine->setDivision(120); // ticks per quarter note
    m_engine->setFileFormat(0); // single track
    m_engine->setTracks(1);
    m_engine->writeToFile(fileName);
}

void SMFBuilder::generate()
{
    QList<float> times = { 0.0,  1.0,  2.0,  3.0,  4.0,  5.0,  6.0,  9.0,  9.0,  9.0, 10.0, 10.0, 10.0, 11.0, 11.0};
    QList<float> notes = {57.0, 59.0, 60.0, 62.0, 64.0, 65.0, 67.0, 57.0, 60.0, 64.0, 65.0, 62.0, 59.0, 64.0, 67.0};
    long quarter_length = 120; // quarter note duration in ticks
    int velocityOn = 100; // forte

    for ( int i=0; i<times.length(); ++i) {
        long event_time = static_cast<long>(times[i] * quarter_length);
        int midi_note = static_cast<int>(notes[i]);
        
        // insert note on event
        MidiEvent noteOn;
        noteOn.absTime = event_time;
        noteOn.status = note_on;
        noteOn.data1 = midi_note;
        noteOn.data2 = velocityOn;
        m_sequence.append(noteOn);

        // insert note off event
        MidiEvent noteOff;
        noteOff.absTime = event_time + quarter_length;
        noteOff.status = note_off;
        noteOff.data1 = midi_note;
        noteOff.data2 = 0; // note off
        m_sequence.append(noteOff);
    }
    m_sequence.sort();
}