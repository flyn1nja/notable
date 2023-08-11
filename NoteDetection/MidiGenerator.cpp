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

constexpr float QUARTER_LENGTH_FACTOR = 1/2.4f; // Quarter note duration factor in ticks
constexpr uint VELOCITY_ON_FACTOR = 127;        // Velocity (volume) factor ratio
constexpr double MIN_AMP_TO_PLAY = 0.2;         // Minimum amplitude for a note to start playing
constexpr double MIN_AMP_TO_STOP = 0.0375;      // Minimum amplitude for a note to continue playing
constexpr int MIN_SAMPLES_COUNT = 4;            // Minimum sample count for a note to be valid



std::string NoteOnToString(const MidiGenerator::MidiEvent& evt);
std::string NoteOffToString(const MidiGenerator::MidiEvent& evt);
std::string DelayToString(int delay, int& len);


void insertNoteOffEvt(MidiGenerator::MidiSequence& sequence, int n, int midiPitch);
void insertNoteOnEvt(MidiGenerator::MidiSequence& sequence, int n, int midiPitch, double amplitude);

static float msPerBeat = 0;

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


// Returns n as a time value in milliseconds.
inline int nAsMS(const int n, const float msPerBeat)
{
    return static_cast<int>(n * msPerBeat + 0.5f); // rounded
}

// Returns n as a time value in seconds.
inline float nAsS(const int n, const float msPerBeat)
{
    return n * msPerBeat / 1000;
}

inline void sort(MidiGenerator::MidiSequence& seq)
{
    std::stable_sort(seq.begin(), seq.end());
}

void MidiGenerator::outputToMIDI(const std::vector<double>& maxs, const std::vector<double>& fmaxs, float duration)
{
    std::cout << std::endl;
    std::cout << "Writing MIDI File!" << std::endl;
    std::cout << std::endl;

    // UNUSED
    const uint UNITS_PER_BEAT = swap(10); // quarter note duration in ticks

    msPerBeat = duration / maxs.size() * 1000;

    std::ofstream midiFile("notable_output.mid", std::ios::out | std::ios::binary);

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
    midiFile.write("\x04\x02\x24\x08\0", 5);
    // midiFile.write(reinterpret_cast<const char *>(&TIME_SIG_PARAMS), 5);
    
    // Tempo
    // const uint TEMPO_PARAMS = swap(0x0927c000);
    midiFile.write("\xFF\x51\x03", 3);
    // midiFile.write("\x09\x27\xC0\0", 4);
    midiFile.write("\x12\x4F\x80\0", 4);
    // midiFile.write(reinterpret_cast<const char *>(&TEMPO_PARAMS), 4);

    midiFile.write("\xFF\x2F\x00", 3); // End of track

    // --> Track 2
    midiFile.write("MTrk", 4);

    float centsOffset = 0.0;
    int prevN = 0;
    int inputsBetween = 0;
    int midiPitch = Pitch::getPitchForFrequency(fmaxs[0], &centsOffset);
    double noteMaxAmp = maxs[0];
    bool isNextSameFreq;
    bool isPlayingNote = false;

    bool hasToPlay = false;
    bool hasToStop = false;
    bool willPlay = false;
    bool willStop = false;

    // TODO : Make a new note when the amplitude doubles for the same frequency (strumming)

    bool skipUpdate = false;

    MidiGenerator::MidiSequence sequence = MidiGenerator::MidiSequence();

    for (int n = 0; n < maxs.size()-1; ++n)
    {
        // if (abs(maxs[n] - noteMaxAmp) <= 0.1 
        //  && maxs[n] < MIN_AMP_TO_PLAY
        //  && n < maxs.size()-MIN_SAMPLES_COUNT  // End of file reached
        //  && !isPlayingNote)
        // {
        //     prevN = n;
        //     continue;
        // }

        isNextSameFreq = abs(midiPitch - Pitch::getPitchForFrequency(fmaxs[n+1])) == 0;

        // Ignore false note changes (cancel brief jumps, go see further)
        if (!isNextSameFreq && isPlayingNote)
        {
            for (int i = 2; i < MIN_SAMPLES_COUNT && i + n < maxs.size()-1; ++i)
            {
                // Ignore if it's a silence / note break
                if (maxs[n+i] < MIN_AMP_TO_STOP) break;

                // If the pitch goes back to original, cancel the note change
                if (Pitch::getPitchForFrequency(fmaxs[n+i]) == midiPitch)
                {
                    isNextSameFreq = true;
                    // skipUpdate = true;
                    break;
                }
            }

            // if (skipUpdate)
            // {
            //     skipUpdate = false;
            //     continue;
            // }
        }
        else if (!willStop)
        {
            // Save the new max amplitude
            if (maxs[n] > noteMaxAmp)
                noteMaxAmp = maxs[n];

            ++inputsBetween;
        }

        hasToPlay = !willPlay                                // Not already playing
                && noteMaxAmp >= MIN_AMP_TO_PLAY &&          // Of signigficant amplitude
                (sequence.size() == 0                        // First note 
                || !isPlayingNote                            // Has reached high enough amplitude
                || !isNextSameFreq);                         // New note (another frequency)

        hasToStop = isPlayingNote &&                             // Has to be playing
                // && n - prevN > MIN_SAMPLES_COUNT &&          // Is long enough in duration
                (maxs[n] < MIN_AMP_TO_STOP                   // Note no longer of high enough amplitude
                || n >= maxs.size()-MIN_SAMPLES_COUNT        // End of file reached
                || !isNextSameFreq);                         // New note (another frequency)


        if (hasToStop && !willStop)
        {
            std::cout << "HasToStop at " << nAsS(n, msPerBeat);

            if (maxs[n] < MIN_AMP_TO_STOP )
                std::cout << " : Note below amplitude threshold ("
                    << maxs[n] << " < " << MIN_AMP_TO_STOP << ")";
            else if (n >= maxs.size()-MIN_SAMPLES_COUNT)
                std::cout << " : End of file reached";
            else if (!isNextSameFreq)
                std::cout << " : Changing freq from " << midiPitch 
                    << " to " << Pitch::getPitchForFrequency(fmaxs[n+1]);
                          
            std::cout << std::endl;
            willStop = true;
            isPlayingNote = false;
        }
        if (hasToPlay && !willPlay)
        {
            std::cout << "HasToPlay at " << nAsS(n, msPerBeat);
            
            if (sequence.size() == 0)
                std::cout << " : First note";
            else if (!isPlayingNote)
                std::cout << " : Note loud enough to play ("
                    << maxs[n] << " > " << MIN_AMP_TO_PLAY << ")";
            else if (!isNextSameFreq)
                std::cout << " : Changing freq from " << midiPitch 
                    << " to " << Pitch::getPitchForFrequency(fmaxs[n+1]);
                     
            std::cout << std::endl;
            willPlay = true;
            isPlayingNote = true;
            midiPitch = Pitch::getPitchForFrequency(fmaxs[n], &centsOffset);

            // Save where the note actually starts playing
            if (sequence.size() == 0)
                prevN = n;
        }

        if (willPlay && willStop)
        {
            if (n-prevN > MIN_SAMPLES_COUNT)
            {
                insertNoteOnEvt(sequence, prevN, midiPitch, noteMaxAmp);
                insertNoteOffEvt(sequence, n, midiPitch);

                std::cout << "Note On at " << nAsS(prevN, msPerBeat) << "  --  Note Off at " 
                << nAsS(n, msPerBeat) << " - Midi: " << midiPitch << std::endl;

                willPlay = false;
                willStop = false;

                prevN = n;
                noteMaxAmp = maxs[n];
            }
            else
            {
                std::cout << "Note at " << nAsS(prevN, msPerBeat) << " - " 
                << nAsS(n, msPerBeat) << " was skipped" << " - Midi: " << midiPitch << std::endl;
                // midiPitch = Pitch::getPitchForFrequency(fmaxs[n], &centsOffset);
                willPlay = false;
                willStop = false;
            }

            std::cout << std::endl;

            inputsBetween = 0;
        }
    }

    //Place the last note (if necessary)
    if (willPlay && willStop)
    {
        insertNoteOnEvt(sequence, prevN, midiPitch, noteMaxAmp);
        insertNoteOffEvt(sequence, maxs.size()-1, midiPitch);

        std::cout << "Note On at " << nAsS(prevN, msPerBeat) << "  --  Note Off at "
        << nAsS(maxs.size()-1, msPerBeat) << " - Midi: " << midiPitch << std::endl;
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
            ++totalLen;
        }

        totalLen += 3;

        if (n + 1 < sequence.size())
        {
            int len = 0;
            outputNotes.append(DelayToString(sequence[n+1].absTime - sequence[n].absTime, len));
            totalLen += len;
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

    int lenCStr = swap(static_cast<int>(outputNotes.length()) + 3);

    // Length of track
    midiFile.write(reinterpret_cast<const char *>(&lenCStr), 4);
    midiFile.write("\00", 1);

    // Track name
    midiFile.write("\xFF\x03", 2);
    midiFile.write("\00\00", 2);

    midiFile.write(outputNotes.c_str(), totalLen);

    std::cout << "End of track" << std::endl;
    midiFile.write("\xFF\x2F\x00", 3); // End of track

    midiFile.close();
}

void MidiGenerator::outputToMIDIPolyphonic(const std::vector<std::vector<double>>& maxs, const std::vector<std::vector<double>>& fmaxs)
{
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


    MidiGenerator::MidiSequence sequence = MidiGenerator::MidiSequence();

    std::vector<bool> isNextSameFreq = std::vector<bool>();
    std::vector<bool> isPlayingNote = std::vector<bool>(12, false);
    std::vector<bool> hasToPlay = std::vector<bool>(12, false);
    std::vector<bool> hasToStop = std::vector<bool>(12, false);
    std::vector<bool> willPlay = std::vector<bool>(12, false);
    std::vector<bool> willStop = std::vector<bool>(12, false);
    std::vector<int> lastN = std::vector<int>(12, 0);
    std::vector<int> midiPitch = std::vector<int>(12, 0);
    std::vector<double> noteMaxAmp = std::vector<double>(12, 0);

    for (int m = 0; m < maxs[0].size()-1; ++m)
    {
        midiPitch[m] = Pitch::getPitchForFrequency(fmaxs[0][m], &centsOffset);
        isNextSameFreq.push_back(fmaxs[1].size()-1 < m 
                                || abs(midiPitch[m] - Pitch::getPitchForFrequency(fmaxs[1][m])) == 0);

        noteMaxAmp[m] = maxs[0][m];
    }

    for (int n = 0; n < maxs.size()-1; ++n)
    {
        std::vector<bool> isPlayingNote(false, maxs[n].size());
        std::vector<bool> hasToPlay(false, maxs[n].size());
        std::vector<bool> hasToStop(false, maxs[n].size());
        std::vector<bool> willPlay(false, maxs[n].size());
        std::vector<bool> willStop(false, maxs[n].size());
        std::vector<int> midiPitch;
        std::vector<double> noteMaxAmp;

        for (int m = 0; m < maxs[n].size(); ++m)
        {
            midiPitch.push_back(Pitch::getPitchForFrequency(fmaxs[n][m], &centsOffset));
            isNextSameFreq.push_back(fmaxs[n+1].size()-1 < m 
                                  || abs(midiPitch[m] - Pitch::getPitchForFrequency(fmaxs[n+1][m])) == 0);
            noteMaxAmp.push_back(maxs[n][m]);
        }

        for (int m = 0; m < maxs[n].size()-1; ++m)
        {

            hasToPlay[m] = !willPlay[m]                                 // Not already playing
                        && noteMaxAmp[m] >= 1 &&                           // Of signigficant amplitude
                        (sequence.size() == 0                           // First note 
                        || !isNextSameFreq[m]                           // New note (another frequency)
                        || !isPlayingNote[m]);                            // Has reached high enough amplitude

            hasToStop[m] = isPlayingNote[m] &&                           // Has to be playing
                        (noteMaxAmp[m] < 0.5                           // Note no longer of high enough amplitude
                        || n+1 == maxs.size()-1                        // End of file reached
                        || !isNextSameFreq[m]                          // New note (another frequency)
                        || maxs[n+1].size()-1 < m);                    // Not enough note on the next  


            if (willPlay[m] && willStop[m])
            {
                insertNoteOnEvt(sequence, lastN[m], midiPitch[m], noteMaxAmp[m]);
                insertNoteOffEvt(sequence, n, midiPitch[m]);

                midiPitch[m] = Pitch::getPitchForFrequency(fmaxs[n+1][m], &centsOffset);

                lastN[m] = n+1;
                noteMaxAmp = maxs[n+1];
                willPlay[m] = false;
                willStop[m] = false;
            }
            else if (willPlay[m] && willStop[m])
            {
                // Test phase to check if it doesn't do some funky stuff (ignores small notes)
                // if (hasToPlay || hasToStop)
                // {
                //     willPlay = false;
                //     willStop = false;
                // }
            }

            if (willStop[m] && !willPlay[m])
            {
                insertNoteOffEvt(sequence, n, midiPitch[m]);
                isPlayingNote[m] = false;
                noteMaxAmp[m] = 0;
                willStop[m] = false;
            }

            if (!willStop[m])
            {
                if (noteMaxAmp < maxs[n])
                {
                    noteMaxAmp = maxs[n];
                    lastN[m] = n+1;
                }
            }

            if (hasToPlay[m])
            {
                // std::cout << "HasToPlay at " << n << std::endl;
                willPlay[m] = true;
                isPlayingNote[m] = true;
                lastN[m] = n+1;
            }
            if (hasToStop[m])
            {
                // std::cout << "HasToStop at " << n << std::endl;
                willStop[m] = true;
                isPlayingNote[m] = false;
            }
        }

    }

    for (int m = 0; m < maxs[maxs.size()-1].size(); ++m)
    {
        //Place the last note (if necessary)
        if (willPlay[m] && willStop[m])
        {
            insertNoteOnEvt(sequence, lastN[m], midiPitch[m], noteMaxAmp[m]);
            insertNoteOffEvt(sequence, maxs.size()-1, midiPitch[m]);
        }
    }

    // Sort the sequence
    sort(sequence);

    using namespace std::string_literals;

    std::string outputNotes = ""s;

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
            int len = 0;
            outputNotes.append(DelayToString(sequence[n+1].absTime - sequence[n].absTime, len));
            totalLen += len;
        }

    }

    int lenCStr = swap(static_cast<int>(outputNotes.length()) + 3);

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
    noteOff.absTime = nAsMS(n, QUARTER_LENGTH_FACTOR * msPerBeat);
    noteOff.status = NOTE_OFF;
    noteOff.data1 = midiPitch;
    noteOff.data2 = 0; // Note off velocity, use 0
    sequence.push_back(noteOff);
}

void insertNoteOnEvt(MidiGenerator::MidiSequence& sequence, int n, int midiPitch, double amplitude)
{
    // Insert note on event
    MidiGenerator::MidiEvent noteOn;
    noteOn.absTime = nAsMS(n, QUARTER_LENGTH_FACTOR * msPerBeat);
    noteOn.status = NOTE_ON;
    noteOn.data1 = midiPitch;
    noteOn.data2 = std::min(static_cast<int>(VELOCITY_ON_FACTOR * amplitude + 0.5), 127);
    sequence.push_back(noteOn);
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

std::string DelayToString(int delay, int& len)
{
    using namespace std::string_literals;
    constexpr int DELAY_INIT = 0x80;

    std::string out = "\x00"s;

    len = 1;

    if (delay >= DELAY_INIT)
    {
        // For each set of 0x80 (128) units delay,
        // increment the deltaTime prefix
        int deltaTimePrefix = DELAY_INIT;

        do
        {
            delay -= DELAY_INIT;
            ++deltaTimePrefix;
        } while (delay >= DELAY_INIT);
        
        out = reinterpret_cast<char *>(&deltaTimePrefix);
        ++len;
        out.append(reinterpret_cast<char *>(&delay), 1);
    }
    else if (delay != 0)
        out = reinterpret_cast<char *>(&delay);

    return out;
}