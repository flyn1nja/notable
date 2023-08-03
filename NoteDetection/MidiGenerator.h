#include <vector>

namespace MidiGenerator
{
    struct MidiEvent
    {
        long absTime;
        int  status;
        int  data1;
        int  data2;

    };

    typedef std::vector<MidiEvent> MidiSequence;

    void outputToMIDI(const std::vector<double>& maxs, const std::vector<double>& fmaxs);

    inline bool operator<(MidiGenerator::MidiEvent & lhs, MidiGenerator::MidiEvent & rhs) { return lhs.absTime < rhs.absTime; }
    inline bool operator<(const MidiGenerator::MidiEvent & lhs, MidiGenerator::MidiEvent & rhs) { return lhs.absTime < rhs.absTime; }
    inline bool operator<(MidiGenerator::MidiEvent & lhs, const MidiGenerator::MidiEvent & rhs) { return lhs.absTime < rhs.absTime; }
    inline bool operator<(const MidiGenerator::MidiEvent & lhs, const MidiGenerator::MidiEvent & rhs) { return lhs.absTime < rhs.absTime; }
}