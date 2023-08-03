#include "MidiGenerator.h"

#include "../constant-q-cpp-master/src/Pitch.h"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdint>
#include <iomanip>

#if BYTE_ORDER == BIG_ENDIAN
#define IS_CORRECT_ENDIAN true
#else
#define IS_CORRECT_ENDIAN false
#endif

constexpr short NOTE_ON = 0b1001;
constexpr short NOTE_OFF = 0b1000;

// const uint UNITS_PER_BEAT = swap(40); // quarter note duration in ticks
constexpr ulong QUARTER_LENGTH_FACTOR = 5; // quarter note duration factor in ticks
constexpr uint VELOCITY_ON_FACTOR = 20;


std::string NoteOnToString(const MidiGenerator::MidiEvent& evt);
std::string NoteOffToString(const MidiGenerator::MidiEvent& evt);
std::string DelayToString(int delay);


void insertNoteOffEvt(MidiGenerator::MidiSequence& sequence, int n, int midiPitch);
void insertNoteOnEvt(MidiGenerator::MidiSequence& sequence, int n, int midiPitch, double amplitude);

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



inline void sort(MidiGenerator::MidiSequence& seq)
{
    std::stable_sort(seq.begin(), seq.end());
}

void MidiGenerator::outputToMIDI(const std::vector<double>& maxs, const std::vector<double>& fmaxs)
{
    // UNUSED
    const uint UNITS_PER_BEAT = swap(10); // quarter note duration in ticks

    std::ofstream midiFile("notable_output.mid", std::ios::out | std::ios::binary );

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
    // midiFile.write(reinterpret_cast<const char *>(&UNITS_PER_BEAT), 2); // unit of time for delta timing. It represents the units per beat. For example, +96 would mean 96 ticks per beat
    midiFile.write("\1\xE0", 2); // number of tracks that follow 
    

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
    // const u_long TIME_SIG_PARAMS = swap(0x0402240800);
    midiFile.write("\xFF\x58\x04", 3);
    midiFile.write("\x04\x02\x18\x08\0", 5);
    // midiFile.write(reinterpret_cast<const char *>(&TIME_SIG_PARAMS), 5);
    
    // Tempo
    // const uint TEMPO_PARAMS = swap(0x0927c000);
    midiFile.write("\xFF\x51\x03", 3);
    midiFile.write("\x09\x27\xC0\0", 4);
    // midiFile.write(reinterpret_cast<const char *>(&TEMPO_PARAMS), 4);

    midiFile.write("\xFF\x2F\x00", 3); // End of track
    
    std::cout << "Track!" << std::endl;

    // --> Track 2
    midiFile.write("MTrk", 4);

    float centsOffset = 0.0;
    int lastN = 0;
    int inputsBetween = 2;
    int midiPitch = Pitch::getPitchForFrequency(fmaxs[0], &centsOffset);
    double noteMaxAmp = maxs[0];
    bool isNextSameFreq;
    bool notePlaying = false;

    bool hasToPlay = false;
    bool hasToStop = false;
    bool willPlay = false;
    bool willStop = false;

    MidiGenerator::MidiSequence sequence = MidiGenerator::MidiSequence();



    for (int n = 0; n < maxs.size()-1; ++n)
    {
        isNextSameFreq = abs(midiPitch - Pitch::getPitchForFrequency(fmaxs[n+1])) == 0;
                        //  && sequence.size() != 0;

        // std::cout << abs(Pitch::getPitchForFrequency(fmaxs[n])
        //              - Pitch::getPitchForFrequency(fmaxs[n+1])) << std::endl;

         


        // if (n+1 == maxs.size()-1 && notePlaying) // End of file reached
        // {
            
        //     insertNoteOffEvt(sequence, n, midiPitch);
        //     notePlaying = false;
        //     break;
        // }
        // else if (isNextSameFreq) // Same note is being played
        // {
        //     if (maxs[n] < 1)
        //     {
        //         if (!notePlaying)
        //         {
        //             ++inputsBetween;
        //             continue;
        //         }
        //         else
        //         {
        //             insertNoteOffEvt(sequence, n, midiPitch);
        //             notePlaying = false;
        //         }
        //     }
        //     else
        //     {
        //         if (maxs[n] > noteMaxAmp)
        //             noteMaxAmp = maxs[n];

        //         ++inputsBetween;
        //     }
        // }
        // else if (!isNextSameFreq          // Changing frequency
        //         && noteMaxAmp >= 1        // Amp is considerably high
        //         && inputsBetween >= 1     // At least 2 samples apart
        //         || !notePlaying && noteMaxAmp >= 1)
        // {
        //     if (notePlaying)
        //         insertNoteOffEvt(sequence, n, midiPitch);
        //     else
        //         notePlaying = true;

        //     // Insert note on event
        //     insertNoteOnEvt(sequence, lastN, midiPitch, std::max(noteMaxAmp, maxs[n]));

        //     // TODO  Mettre le cents offset dans le pitch?
        //     midiPitch = Pitch::getPitchForFrequency(fmaxs[n], &centsOffset);

        //     lastN = n;
        //     noteMaxAmp = maxs[n];
        //     inputsBetween = 0;
        // }
        // else 
        // {
        //     ++inputsBetween;
        // }


       hasToPlay = !willPlay                                 // Not already playing
                && noteMaxAmp >= 1 &&                        // Of signigficant amplitude
                (sequence.size() == 0                        // First note 
                || !isNextSameFreq                           // New note (another frequency)
                || !notePlaying);                            // Has reached high enough amplitude

       hasToStop = notePlaying &&              // Has to be playing
                   (noteMaxAmp < 0.5                         // Note no longer of high enough amplitude
                 || n+1 == maxs.size()-1                     // End of file reached
                 || !isNextSameFreq);                        // New note (another frequency)


        if (willPlay && willStop)
        //  && inputsBetween >= 2 && hasToPlay)
        {
            insertNoteOnEvt(sequence, lastN, midiPitch, noteMaxAmp);
            insertNoteOffEvt(sequence, n, midiPitch);

            // notePlaying = true;

            std::cout << "Note On at " << lastN << "  --  Note Off at " << n << std::endl;

            midiPitch = Pitch::getPitchForFrequency(fmaxs[n+1], &centsOffset);

            lastN = n+1;
            noteMaxAmp = maxs[n+1];
            inputsBetween = 0;
            willPlay = false;
            willStop = false;
        }
        else if (willPlay && willStop)
        {
            // Test phase to check if it doesn't do some funky stuff (ignores small notes)
            ++inputsBetween;
            // if (hasToPlay || hasToStop)
            // {
            //     willPlay = false;
            //     willStop = false;
            // }
        }

        if (willStop && !willPlay)
        {
            insertNoteOffEvt(sequence, n, midiPitch);
            notePlaying = false;
            noteMaxAmp = 0;
            willStop = false;
        }

        if (!willStop)
        {
            if (noteMaxAmp < maxs[n])
            {
                noteMaxAmp = maxs[n];
                lastN = n;
            }
            ++inputsBetween;
        }

        if (hasToPlay)
        {
            std::cout << "HasToPlay at " << n << std::endl;
            willPlay = true;
            notePlaying = true;
            lastN = n;
        }
        if (hasToStop)
        {
            std::cout << "HasToStop at " << n << std::endl;
            willStop = true;
            notePlaying = false;
        }
    }

    //Place the last note (if necessary)
    if (willPlay && willStop)
    {
        insertNoteOnEvt(sequence, lastN, midiPitch, noteMaxAmp);
        insertNoteOffEvt(sequence, maxs.size()-1, midiPitch);

        // notePlaying = true;

        std::cout << "Note On at " << lastN << "  --  Note Off at " << maxs.size()-1 << std::endl;
    }

    // Sort the sequence
    sort(sequence);

    using namespace std::string_literals;

    std::string outputNotes = ""s;

    // std::vector<int> values = std::vector<int>();
#if 1
    int totalLen = 0;

    for (int n = 0; n < sequence.size(); ++n)
    {

        if (sequence[n].status == 9)
            outputNotes.append(NoteOnToString(sequence[n]));
        else
        {
            outputNotes.append(NoteOffToString(sequence[n]));
            outputNotes.append("\0"s);
            ++totalLen;
        }

        totalLen += 3;


        if (n + 1 < sequence.size() && sequence[n+1].status == 8)
        {
            outputNotes.append(DelayToString(sequence[n+1].absTime - sequence[n].absTime));
            totalLen += 2;
        }

    }
#else
    int timeDelta = 0;
    int outputVal = 0x0;
    int delay = 0x0;
    char out[11];
    std::vector<bool> hasDelay = std::vector<bool>();
    bool followedByDelay = false;

    for (int n = 0; n < sequence.size(); ++n)
    {

        if (n + 2 < sequence.size() && sequence[n+2].status == 8)
            timeDelta = sequence[n+1].absTime - sequence[n].absTime;
        else
            timeDelta = 0;


        outputVal = sequence[n].status << 0x14;
        outputVal |= sequence[n].data1 << 0x8;
        outputVal |= sequence[n].data2;
        
        followedByDelay = timeDelta > 0;
        hasDelay.push_back(followedByDelay);

        if (followedByDelay)
        {
            delay = 0x8100;
            delay |= std::min(timeDelta, 0xFF);
        }
        else
            delay = 0;

        std::printf("val: %06lx  delay: %04lx\t", outputVal, delay);

        // outputVal = swap(outputVal);
        // delay = swap(delay);


        // track_event = <v_time> + <midi_event> | <meta_event> | <sysex_event>
        int d = std::snprintf(out, 7, "%06lx", outputVal);
        outputNotes.append(out);

        if (timeDelta > 0)
        {
            std::snprintf(out+d, 5, "%04lx", delay);
            outputNotes.append(out);
        }
    
        outputVal = 0;
        delay = 0;
    }

    hasDelay.push_back(false);

    int lenCStr = swap(static_cast<int>(outputNotes.length()));

    // Length of track
    midiFile.write(reinterpret_cast<const char *>(&lenCStr), 4);

    // Track name
    midiFile.write("\xFF\x03", 2);
    midiFile.write("\00\00", 2);

    std::cout << "Writing notes: " << outputNotes << std::endl;

    std::stringstream ss;
    int chunkSize = 3;
    int currSize = 0;
    int val = 0;
    int pos = 0;

    for (int c = 0; c < outputNotes.length(); c+=2)
    {
        ss << outputNotes.substr(c,2);
        ++currSize;
        if (c > 0 && currSize == chunkSize)
        {
            ss >> std::hex >> val;
            std::cout << std::hex << val << std::endl;
            val = swap(val);
            midiFile.write(reinterpret_cast<const char *>(&val), chunkSize);

            ++pos;
            chunkSize = hasDelay[pos] ? 2 : 3;
            ss.clear();
            val = 0;
            currSize = 0;
        }
        else if (outputNotes.length() - c < chunkSize)
            ss >> val;
        
    }

    val = swap(val);

    midiFile.write(reinterpret_cast<const char *>(&val), chunkSize);
#endif
    // ss << std::hex << outputNotes;
    // midiFile.write(, lenCStr);

    int lenCStr = swap(static_cast<int>(outputNotes.length()));

    // Length of track
    midiFile.write(reinterpret_cast<const char *>(&lenCStr), 4);
    midiFile.write("\00", 1);

    // Track name
    midiFile.write("\xFF\x03", 2);
    midiFile.write("\00\00", 2);

    midiFile.write(outputNotes.c_str(), totalLen); // End of track

    std::cout << "End of track" << std::endl;
    midiFile.write("\xFF\x2F\x00", 3); // End of track

    midiFile.close();
}

void insertNoteOffEvt(MidiGenerator::MidiSequence& sequence, int n, int midiPitch)
{
    // Insert note off event
    MidiGenerator::MidiEvent noteOff;
    noteOff.absTime = (n * QUARTER_LENGTH_FACTOR)/50;
    noteOff.status = NOTE_OFF;
    noteOff.data1 = midiPitch;
    noteOff.data2 = 0; // Note off velocity, use 0
    sequence.push_back(noteOff);
}

void insertNoteOnEvt(MidiGenerator::MidiSequence& sequence, int n, int midiPitch, double amplitude)
{
    // Insert note on event
    MidiGenerator::MidiEvent noteOn;
    noteOn.absTime = (n * QUARTER_LENGTH_FACTOR)/50;
    noteOn.status = NOTE_ON;
    noteOn.data1 = midiPitch;
    noteOn.data2 = std::min(VELOCITY_ON_FACTOR * amplitude, 127.0);
    sequence.push_back(noteOn);
}

template<typename T>
std::string intToHex(T integer, int w=sizeof(T)*2)
{
  std::stringstream ss;
  ss << "0x" << std::setfill ('\0') << std::setw(w) 
         << std::hex << integer;
  return ss.str();
}

std::string NoteOnToString(const MidiGenerator::MidiEvent& evt)
{
    using namespace std::string_literals;

    int pitch = evt.data1;
    int velocity = evt.data2;
    std::string out = "\x90"s;
    out.append(reinterpret_cast<char *>(&pitch),1);
    out.append(reinterpret_cast<char *>(&velocity),1);

    return out;
    // return out.append(reinterpret_cast<const char *>(&pitch)).append(reinterpret_cast<const char *>(&velocity));
}

std::string NoteOffToString(const MidiGenerator::MidiEvent& evt)
{
    using namespace std::string_literals;

    int pitch = evt.data1;
    std::string out = "\x80"s;
    out.append(reinterpret_cast<char *>(&pitch));
    out.append("\0"s);

    return out;
}

std::string DelayToString(int delay)
{
    using namespace std::string_literals;

    std::string out = "\x81"s;
    delay = std::min(delay, 0xFF);
    out.append(reinterpret_cast<char *>(&delay),1);

    return out;
}