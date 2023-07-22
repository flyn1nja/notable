#include "MidiGenerator.h"

#include "../constant-q-cpp-master/src/Pitch.h"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdint>

#if BYTE_ORDER == BIG_ENDIAN
#define IS_CORRECT_ENDIAN true
#else
#define IS_CORRECT_ENDIAN false
#endif

constexpr short NOTE_ON = 0b1001;
constexpr short NOTE_OFF = 0b1000;

inline void sort(MidiGenerator::MidiSequence& seq)
{
    std::stable_sort(seq.begin(), seq.end());
}

// Swaps the bytes if needed
template <typename T> 
T swap(const T& arg)
{
    #if (IS_CORRECT_ENDIAN)
    return arg;             // no byte-swapping needed

    #else                     // swap bytes
    T ret;

    char* dst = reinterpret_cast<char*>(&ret);
    const char* src = reinterpret_cast<const char*>(&arg + 1);

    for (size_t i = 0; i < sizeof(T); i++)
        *dst++ = *--src;

    return ret;
    #endif
}

// void SMFBuilder::trackHandler(int )
// {
//     // meta events
//     m_engine->writeBpmTempo(0, 100); // m.m.=100 (andante)
//     m_engine->writeTimeSignature(0, 4, 2, 18, 8);  // ts = 4/4
//     m_engine->writeKeySignature(0, 0, major_mode); // C major
//     m_engine->writeMidiEvent(0, program_chng, 0, 0);  // grand piano
//     // note events
//     long last_time = 0;
//     for(auto ev : m_sequence)
//     {
//         long delta = ev.absTime - last_time;
//         last_time = ev.absTime;
//         m_engine->writeMidiEvent(delta,  ev.status, 0, ev.data1, ev.data2);
//     }
//     // final event
//     m_engine->writeMetaEvent(0, end_of_track); 
// }

// void SMFBuilder::run(QString fileName)
// {
//     m_engine->setDivision(120); // ticks per quarter note
//     m_engine->setFileFormat(0); // single track
//     m_engine->setTracks(1);
//     m_engine->writeToFile(fileName);
// }

// void MidiGenerator::generate()
// {
//     QList<float> times = { 0.0,  1.0,  2.0,  3.0,  4.0,  5.0,  6.0,  9.0,  9.0,  9.0, 10.0, 10.0, 10.0, 11.0, 11.0};
//     QList<float> notes = {57.0, 59.0, 60.0, 62.0, 64.0, 65.0, 67.0, 57.0, 60.0, 64.0, 65.0, 62.0, 59.0, 64.0, 67.0};
//     long quarter_length = 120; // quarter note duration in ticks
//     int VELOCITY_ON_FACTOR = 100; // forte

//     for ( int i=0; i<times.length(); ++i) {
//         long event_time = static_cast<long>(times[i] * quarter_length);
//         int midi_note = static_cast<int>(notes[i]);
        
//         // insert note on event
//         MidiGenerator::MidiEvent noteOn;
//         noteOn.absTime = event_time;
//         noteOn.status = note_on;
//         noteOn.data1 = midi_note;
//         noteOn.data2 = velocityOn;
//         m_sequence.append(noteOn);

//         // insert note off event
//         MidiGenerator::MidiEvent noteOff;
//         noteOff.absTime = event_time + quarter_length;
//         noteOff.status = note_off;
//         noteOff.data1 = midi_note;
//         noteOff.data2 = 0; // note off
//         m_sequence.append(noteOff);
//     }
//     m_sequence.sort();
// }

void MidiGenerator::outputToMIDI(const std::vector<double>& maxs, const std::vector<double>& fmaxs)
{
    const ulong UNITS_PER_BEAT = swap(40); // quarter note duration in ticks
    constexpr ulong QUARTER_LENGTH_FACTOR = 1; // quarter note duration factor in ticks
    constexpr uint VELOCITY_ON_FACTOR = 20;

    std::ofstream midiFile("notable_output.mid", std::ios::out | std::ios::binary );
    float centsOffset = 0.0;
    int midiPitch = 0;
    int lastN = 0;
    int inputsBetween = 2;
    float lastMax = maxs[0];
    float lastF = -1.0f;
    float noteMaxAmp = maxs[0];
    bool notePlaying = false;

    MidiGenerator::MidiSequence sequence = MidiGenerator::MidiSequence();


    // * Write the header
    // --> Header Chunk
    midiFile.write("MThd", 4); // the literal string MThd, or in hexadecimal notation: 0x4d546864. These four characters at the start of the MIDI file indicate that this is a MIDI file. 
    // midiFile.write(reinterpret_cast<const char *>(u_int8_t(6)), 4); // length of the header chunk (always 6 bytes long -- the size of the next three fields). 
    // midiFile.write(reinterpret_cast<const char *>(u_int8_t(0)), 2); // file format. 0 = single track file format 
    // midiFile.write(reinterpret_cast<const char *>(u_int8_t(1)), 2); // number of tracks that follow 
    // // ! Revoir ! //
    // midiFile.write(reinterpret_cast<const char *>(u_int8_t(96)), 2); // unit of time for delta timing. It represents the units per beat. For example, +96 would mean 96 ticks per beat
    midiFile.write("\0\0\0\6", 4); // length of the header chunk (always 6 bytes long -- the size of the next three fields). 
    midiFile.write("\0\1", 2); // file format. 0 = single track file format, 1 = multi track file format 
    midiFile.write("\0\2", 2); // number of tracks that follow 
    // ! Revoir ! //
    midiFile.write(reinterpret_cast<const char *>(&UNITS_PER_BEAT), 2); // unit of time for delta timing. It represents the units per beat. For example, +96 would mean 96 ticks per beat
    
    // Initial info track
    midiFile.write("MTrk", 4);
    midiFile.write("\0\0\0\x17", 4);

    midiFile.write("\00", 1);

    // // ;)
    // midiFile.write("\xFF\1\x31", 3);
    // midiFile.write("Generated With Notable, By Alexis Giguere-Lebel\0\0", 49);

    // Track name
    midiFile.write("\xFF\03", 2);
    midiFile.write("\0\0", 2);

    // Time signature
    const u_long TIME_SIG_PARAMS = swap(0x0402240800);
    midiFile.write("\xFF\x58\x04", 3);
    midiFile.write(reinterpret_cast<const char *>(&TIME_SIG_PARAMS), 5);
    
    // Tempo
    const uint TEMPO_PARAMS = swap(0x0927c000);
    midiFile.write("\xFF\x51\x03", 3);
    midiFile.write(reinterpret_cast<const char *>(&TEMPO_PARAMS), 4);

    midiFile.write("\xFF\x2F\x00", 3); // End of track
    
    std::cout << "Track!" << std::endl;

    // --> Track 2
    midiFile.write("MTrk", 4);

    bool isNextSameFreq;

    for (int n = 0; n < maxs.size()-1; ++n)
    {
        isNextSameFreq = n != 0 && abs(Pitch::getPitchForFrequency(fmaxs[n])
                         - Pitch::getPitchForFrequency(fmaxs[n+1])) != 0;

        if (maxs[n] >= 1 && abs(maxs[n] - lastMax) > 0.01)
        {
            if (!isNextSameFreq && lastN != n)
            {
                if (notePlaying)
                {
                    // Insert note off event
                    MidiGenerator::MidiEvent noteOff;
                    noteOff.absTime = n * QUARTER_LENGTH_FACTOR;
                    noteOff.status = NOTE_OFF;
                    noteOff.data1 = midiPitch;
                    noteOff.data2 = 0; // Note off velocity, use 0
                    sequence.push_back(noteOff);
                }
                else
                    notePlaying = true;

                // TODO  Mettre le cents offset dans le pitch?
                midiPitch = Pitch::getPitchForFrequency(fmaxs[n], &centsOffset);

                long eventTime = static_cast<long>(n - lastN);

                // Insert note on event
                MidiGenerator::MidiEvent noteOn;
                noteOn.absTime = lastN * QUARTER_LENGTH_FACTOR;
                noteOn.status = NOTE_ON;
                noteOn.data1 = midiPitch;
                noteOn.data2 = VELOCITY_ON_FACTOR * noteMaxAmp;
                sequence.push_back(noteOn);

                lastN = n+1;
                lastF = fmaxs[n];
                noteMaxAmp = 0;
                inputsBetween = 0;

                // midiFile << ;
            }
            else 
            {
                if (notePlaying && maxs[n] < 1 && isNextSameFreq)
                {
                    notePlaying = false;

                    // Insert note off event
                    MidiGenerator::MidiEvent noteOff;
                    noteOff.absTime = n * QUARTER_LENGTH_FACTOR;
                    noteOff.status = NOTE_OFF;
                    noteOff.data1 = midiPitch;
                    noteOff.data2 = 0; // Note off velocity, use 0
                    sequence.push_back(noteOff);
                }
                else
                {
                    if (isNextSameFreq && maxs[n] > noteMaxAmp)
                        noteMaxAmp = maxs[n];
                    else if (!isNextSameFreq)
                    {
                        lastF = fmaxs[n];
                        lastN = n+1;
                    }
                }

                ++inputsBetween;
            }
        }
    }

    // Sort the sequence
    sort(sequence);

    std::string outputNotes;
    long prevTime = 0l;
    int timeDelta = 0;
    ulong outputVal = 0x0;
    ulong delay = 0x0;
    char out[16];

    for (int n = 0; n < sequence.size(); ++n)
    {
        timeDelta = sequence[n].absTime - prevTime;


        outputVal = sequence[n].status << 0x14;
        outputVal |= sequence[n].data1 << 0x8;
        outputVal |= sequence[n].data2;
        
        if (timeDelta > 0)
        {
            delay = 0x8100;
            delay |= timeDelta;
        }
        else
            delay = 0;

        std::printf("val: %06lx  delay: %04lx\t", outputVal, delay);

        // outputVal = swap(outputVal);
        // delay = swap(delay);


        // track_event = <v_time> + <midi_event> | <meta_event> | <sysex_event>
        int d = std::snprintf(out, 6, "%06lx", outputVal);
        outputNotes.append(out);

        if (timeDelta > 0)
        {
            std::snprintf(out+d, 4, "%04lx", delay);
            outputNotes.append(out);
        }
    
        outputVal = 0;
        delay = 0;
        prevTime = sequence[n].absTime;
    }

    ulong lenCStr = swap(outputNotes.length());

    // Length of track
    midiFile.write(reinterpret_cast<const char *>(&lenCStr), 4);

    // Track name
    midiFile.write("\xFF\x03", 2);
    midiFile.write("\00\00", 2);

    std::cout << "Writing notes: " << outputNotes << std::endl;

    std::stringstream ss;
    constexpr int CHUNK_SIZE = 6;
    uint val = 0;
    int remainder = 0;

    for (int c = 0; c < outputNotes.length(); c+=2)
    {
        ss << outputNotes.substr(c,2);
        if (c > 0 && c % CHUNK_SIZE == 0)
        {
            ss >> std::hex >> val;
            std::cout << val << std::endl;
            val = swap(val);
            midiFile.write(reinterpret_cast<const char *>(&val), CHUNK_SIZE);
            val = 0;
            ss.clear();
        }
        else if (outputNotes.length() - c < CHUNK_SIZE)
        {
            ss >> val;
            ++remainder;
        }
        
    }

    val = swap(val);
    midiFile.write(reinterpret_cast<const char *>(&val), remainder);

    // ss << std::hex << outputNotes;
    // midiFile.write(, lenCStr);



    std::cout << "End of track" << std::endl;
    midiFile.write("\xFF\x2F\x00", 3); // End of track

    midiFile.close();
}