/*
 * File stores many constants useful for the creation
 * of midi files
 * reference for midi specification from: https://www.csie.ntu.edu.tw/~r92092/ref/midi/
 * Author: Charlie Street
 */

#ifndef MIDIMAKER_MIDI_H
#define MIDIMAKER_MIDI_H

#include <vector>
#include <fstream>

//some types i need made easier
typedef unsigned char byte;
typedef std::vector<byte> VLQ;
typedef unsigned __int32 uint32_t;
typedef unsigned __int16 uint16_t;

const std::vector<byte> MThd = {'M','T','h','d'}; //two types to identify chunk
const std::vector<byte> MTrk = {'M','T','r','k'};

/*lots of command codes for midi
 *i am assuming the use of midi channel 0
 */
#define F0_SYSEX 0xF0
#define F7_SYSEX 0xF7
#define META_EVENT 0xFF
#define NOTE_OFF 0x80
#define NOTE_ON 0x90
#define POLYPHONIC_KEY_PRESSURE 0xA0
#define CONTROLLER_CHANGE 0xB0
#define PROGRAM_CHANGE 0xC0
#define CHANNEL_KEY_PRESSURE 0xD0
#define PITCH_BEND 0xE0
#define CHANNEL_MODE 0xB0



//class to represent a generic chunk in a midi file
class MidiChunk {
protected:
    std::vector<byte> type; //will be MThd or MTrk
    uint32_t length;
    std::vector<byte> data;

    void addData(std::vector<byte>); //adds data and increases length appropriately
    void addByte(byte); //same as add data but for single byte

public:
    explicit MidiChunk(std::vector<byte>);
    std::vector<byte> getType();
    uint32_t getLength();
    std::vector<byte> getData();

};


//class to represent the first chunk of a midi file, the header
class MidiHeader : public MidiChunk {
public:
    MidiHeader(uint16_t format, uint16_t tracks, uint16_t division);
};

/* class to represent a track chunk of a midi file
 * i.e. where most of the information is stored
 */
class MidiTrack : public MidiChunk {
public:
    MidiTrack();

    //system exclusive events
    void f0Sysex(unsigned int, std::vector<byte>);
    void f7Sysex(unsigned int, std::vector<byte>);

    //meta events
    void sequenceNum(unsigned int, uint16_t);
    void textEvent(unsigned int, std::string);
    void copyrightNotice(unsigned int, std::string);
    void trackName(unsigned int, std::string);
    void instrumentName(unsigned int, std::string);
    void lyric(unsigned int, std::string);
    void marker(unsigned int, std::string);
    void cuePoint(unsigned int, std::string);
    void midiChannelPrefix(unsigned int, byte);
    void endOfTrack(unsigned int);
    void setTempo(unsigned int, uint32_t); //note, actually 24 bit binary
    void smtpeOffset(unsigned int, uint32_t, byte);
    void timeSignature(unsigned int, byte, byte, byte, byte);
    void keySignature(unsigned int, byte, byte);
    void sequencerSpecificMetaEvent(unsigned int, std::vector<byte>,std::vector<byte>);


    //midi channel voice messages
    void noteOff(unsigned int, byte, byte);
    void noteOn(unsigned int, byte, byte);
    void polyphonicKeyPressure(unsigned int, byte, byte);
    void controllerChange(unsigned int, byte, byte);
    void programChange(unsigned int, byte);
    void channelKeyPressure(unsigned int, byte);
    void pitchBend(unsigned int, byte, byte);

    //midi channel mode messages
    void allSoundOff(unsigned int);
    void resetAllControllers(unsigned int);
    void localControl(unsigned int, byte);
    void allNotesOff(unsigned int);
    void omniModeOff(unsigned int);
    void omniModeOn(unsigned int);
    void monoModeOn(unsigned int, byte);
    void polyModeOn(unsigned int);


};

//class to represent a whole midi file
class MidiFile {
private:

    std::ofstream midiFile;
    void writeChunk(MidiChunk);

public:
    MidiFile(std::string, MidiHeader, MidiTrack);
    void addTrack(MidiTrack);
    void closeFile();

};

//--UTILITY FUNCTIONS--

/* function to convert an unsigned integer into
 * a variable length quantity
 * where the first ms n-1 bytes have msb as 1
 * and the ls byte has the lsb as 0
 * byte 0 is most significant
 * as that is the byte we want in the midi file
 */
VLQ uintToVLQ(unsigned int);

//converts a string into a byte vector
std::vector<byte> fromString(std::string);

/* used to extract bytes out of an unsigned integer
 * convention is least significant byte is byte 0
 */
byte getNthByte16(uint16_t, unsigned int);
byte getNthByte32(uint32_t, unsigned int);

#endif //MIDIMAKER_MIDI_H
