/*
*    Title: Midi-Maker
*    Author: Charlie Street "charlie1329"
*    Date: 26.9.2017
*    Code version: commit e06c37903d821a4923cdc38cf025fec889a88ed8
*    Availability: https://github.com/charlie1329/Midi-Maker
*/

/*
 * File implements the ability to
 * writing raw midi files
 * reference for midi specification from: https://www.csie.ntu.edu.tw/~r92092/ref/midi/
 * Author: Charlie Street
 */

#include <iostream>
#include <iterator>
#include <utility>
#include <algorithm>
#include "midi.h"


//--MIDI CHUNK CLASS IMPLEMENTATIONS--

MidiChunk::MidiChunk(std::vector<byte> type) : type(std::move(type)), length(0) {}

/* credit for finding how to concatenate two vectors in c++
 * goes to Joe/Andreas Brink at:
 * https://stackoverflow.com/questions/2551775/appending-a-vector-to-a-vector
 */
void MidiChunk::addData(std::vector<byte> toAdd) {
    data.insert(std::end(data),std::begin(toAdd),std::end(toAdd)); //append to member vector
    length += toAdd.size(); //increase length appropriately
}

void MidiChunk::addByte(byte toAdd) {
    data.push_back(toAdd);
    length++;
}

//get methods
std::vector<byte> MidiChunk::getType() {
    return type;
}

uint32_t MidiChunk::getLength() {
    return length;
}

std::vector<byte> MidiChunk::getData(){
    return data;
}

//--MIDI HEADER CLASS IMPLEMENTATIONS--

/* the header chunk should always have 6 bytes of data, hence everything can just be done in the constructor
 * format is midi file format (should usually be 0)
 * tracks is number of track chunks in midi file
 * division is default unit of delta time for midi file
 * if bit 15 is 0 then the rest is ticks per quarter note
 * if bit 15 is 1 then 14..8 are -frames/second and 7..0 are ticks/frame
 */
MidiHeader::MidiHeader(uint16_t format, uint16_t tracks, uint16_t division) : MidiChunk(MThd) {
    std::vector<byte> headerData;
    headerData.push_back(getNthByte16(format,1));
    headerData.push_back(getNthByte16(format,0));
    headerData.push_back(getNthByte16(tracks,1));
    headerData.push_back(getNthByte16(tracks,0));
    headerData.push_back(getNthByte16(division,1));
    headerData.push_back(getNthByte16(division,0));
    addData(headerData);
}

//--MIDI TRACK CLASS IMPLEMENTATIONS--
//--NOTE: Any function with a delta parameter represents
//--Number of ticks from previous event
//--And is a Variable Length Quantity
//--Though in my functions they are internally converted
//--To VLQs for ease of use

MidiTrack::MidiTrack() : MidiChunk(MTrk) {}

//--Sysex messages--

//FO <sysex_data> is sent as a midi message
void MidiTrack::f0Sysex(unsigned int delta, std::vector<byte> sysexData) {
    addData(uintToVLQ(delta));
    addByte(F0_SYSEX);
    addData(uintToVLQ(sysexData.size()));
    addData(sysexData);
}

//<any_data> is sent as a midi message
void MidiTrack::f7Sysex(unsigned int delta, std::vector<byte> anyData) {
    addData(uintToVLQ(delta));
    addByte(F7_SYSEX);
    addData(uintToVLQ(anyData.size()));
    addData(anyData);
}

//--Meta Events--

//must only occur at the start of a track (optional)
void MidiTrack::sequenceNum(unsigned int delta, uint16_t seq) {
    addData(uintToVLQ(delta));
    addData({META_EVENT, 0x00, 0x02,
             getNthByte16(seq, 1),
             getNthByte16(seq, 0)});
}

//annotates with arbitrary text
void MidiTrack::textEvent(unsigned int delta, std::string str) {
    addData(uintToVLQ(delta));
    addData({META_EVENT, 0x01});
    std::vector<byte> strVec = fromString(std::move(str));
    addData(uintToVLQ(strVec.size()));
    addData(strVec);
}

//adds copyright to file
void MidiTrack::copyrightNotice(unsigned int delta, std::string str) {
    addData(uintToVLQ(delta));
    addData({META_EVENT, 0x02});
    std::vector<byte> strVec = fromString(std::move(str));
    addData(uintToVLQ(strVec.size()));
    addData(strVec);
}

//gives a name to this track chunk
void MidiTrack::trackName(unsigned int delta, std::string str) {
    addData(uintToVLQ(delta));
    addData({META_EVENT, 0x03});
    std::vector<byte> strVec = fromString(std::move(str));
    addData(uintToVLQ(strVec.size()));
    addData(strVec);
}

//adds a name to a particular instrument
void MidiTrack::instrumentName(unsigned int delta, std::string str) {
    addData(uintToVLQ(delta));
    addData({META_EVENT, 0x04});
    std::vector<byte> strVec = fromString(std::move(str));
    addData(uintToVLQ(strVec.size()));
    addData(strVec);
}

//adds lyrics to the track
void MidiTrack::lyric(unsigned int delta, std::string str) {
    addData(uintToVLQ(delta));
    addData({META_EVENT, 0x05});
    std::vector<byte> strVec = fromString(std::move(str));
    addData(uintToVLQ(strVec.size()));
    addData(strVec);
}

//marks a significant point in the sequence
void MidiTrack::marker(unsigned int delta, std::string str) {
    addData(uintToVLQ(delta));
    addData({META_EVENT, 0x06});
    std::vector<byte> strVec = fromString(std::move(str));
    addData(uintToVLQ(strVec.size()));
    addData(strVec);
}

//used to include stage cues etc.
void MidiTrack::cuePoint(unsigned int delta, std::string str) {
    addData(uintToVLQ(delta));
    addData({META_EVENT, 0x07});
    std::vector<byte> strVec = fromString(std::move(str));
    addData(uintToVLQ(strVec.size()));
    addData(strVec);
}

//associate all following meta-events and sysex-events with specified midi channel
void MidiTrack::midiChannelPrefix(unsigned int delta, byte channel) {
    addData(uintToVLQ(delta));
    addData({META_EVENT, 0x20, 0x01, channel});
}

//non optional event to give track a clearly defined length
void MidiTrack::endOfTrack(unsigned int delta) {
    addData(uintToVLQ(delta));
    addData({META_EVENT, 0x2F, 0x00});
}

//sets the tempo in microseconds per quarter note
//if not specified, the default tempo is 120 bpm, which is ttttt = 50000
void MidiTrack::setTempo(unsigned int delta, uint32_t tempo) {
    addData(uintToVLQ(delta));
    addData({META_EVENT, 0x51, 0x03,
             getNthByte32(tempo, 2),
             getNthByte32(tempo, 1),
             getNthByte32(tempo, 0)});
}

//optional event that specifies smtpe time at which the track is to start
void MidiTrack::smtpeOffset(unsigned int delta, uint32_t time, byte frac) {
    addData(uintToVLQ(delta));
    addData({META_EVENT, 0x54, 0x05,
             getNthByte32(time, 3),
             getNthByte32(time, 2),
             getNthByte32(time, 1),
             getNthByte32(time, 0),
             frac});
}

/* Time signature of form num/2^dom
 * cpt is midi clocks per metronome tick
 * bb is number of 1/32 notes per 34 midi clocks (8 is standard)
 */
void MidiTrack::timeSignature(unsigned int delta, byte num, byte dom, byte cpt, byte bb) {
    addData(uintToVLQ(delta));
    addData({META_EVENT, 0x58, 0x04, num, dom, cpt, bb});
}

/*key signature, expressed as number of sharps/flats, and major/minor flag
 * -7 = 7 flats, 0 = middle C, +7 = 7 sharps
 * 0 = major, 1 = minor
 */
void MidiTrack::keySignature(unsigned int delta, byte sf, byte mi) {
    addData(uintToVLQ(delta));
    addData({META_EVENT, 0x59, 0x02, sf, mi});
}

/* midi-file equivalent of system exclusive message
 * id is 1 or 3 bytes representing the manufactures id
 */
void MidiTrack::sequencerSpecificMetaEvent(unsigned int delta, std::vector<byte> id, std::vector<byte> data) {
    addData(uintToVLQ(delta));
    addData({META_EVENT, 0x7F});
    addData(uintToVLQ(id.size() + data.size()));
    addData(id);
    addData(data);
}

//--MIDI CHANNEL VOICE MESSAGES--

/* normally sent when a key is released
 * 3C = Middle C
 * 40 = Defualt Velocity
 */
void MidiTrack::noteOff(unsigned int delta, byte note, byte vel) {
    addData(uintToVLQ(delta));
    addData({NOTE_OFF, note, vel});
}

/* normally sent when a key is pressed
 *
 */
void MidiTrack::noteOn(unsigned int delta, byte note, byte vel) {
    addData(uintToVLQ(delta));
    addData({NOTE_ON, note, vel});
}


//sent when there is a change in pressure being applied to a key
void MidiTrack::polyphonicKeyPressure(unsigned int delta, byte key, byte pressure) {
    addData(uintToVLQ(delta));
    addData({POLYPHONIC_KEY_PRESSURE, key, pressure});
}

/* sent when a change is made in a footswitch,
 * expression pedal, slider, or other controller
 */
void MidiTrack::controllerChange(unsigned int delta, byte num, byte val) {
    addData(uintToVLQ(delta));
    addData({CONTROLLER_CHANGE, num, val});
}

/* used to change the instrument to be played
 * when a note on message is received
 */
void MidiTrack::programChange(unsigned int delta, byte prog) {
    addData(uintToVLQ(delta));
    addData({PROGRAM_CHANGE, prog});
}

/* Also known as aftertouch.
 * This message is sent when there is a change in the
 * overall pressure being applied to the keyboard
 */
void MidiTrack::channelKeyPressure(unsigned int delta, byte pressure) {
    addData(uintToVLQ(delta));
    addData({CHANNEL_KEY_PRESSURE, pressure});
}

//Sent when a change is mafe in a pitch-bender level
void MidiTrack::pitchBend(unsigned int delta, byte lsb, byte msb) {
    addData(uintToVLQ(delta));
    addData({PITCH_BEND, lsb, msb});
}

//--MIDI CHANNEL MODE MESSAGES

//Turn off all sound including envelopes and reverb effects
void MidiTrack::allSoundOff(unsigned int delta) {
    addData(uintToVLQ(delta));
    addData({CHANNEL_MODE, 0x78, 0x00});
}

//reset all controllers to their default positions
void MidiTrack::resetAllControllers(unsigned int delta) {
    addData(uintToVLQ(delta));
    addData({CHANNEL_MODE, 0x79, 0x00});
}

/* disconnect (or reconnect) the keyboard and
 * sound generator in a MIDI synthesizer
 */
void MidiTrack::localControl(unsigned int delta, byte connect) {
    addData(uintToVLQ(delta));
    addData({CHANNEL_MODE, 0x7A, connect});
}

/* turns off all notes for which a
 * note on message has been received
 */
void MidiTrack::allNotesOff(unsigned int delta) {
    addData(uintToVLQ(delta));
    addData({CHANNEL_MODE, 0x7B, 0x00});
}

//receiver should respong to voice messages on it's basic channel
void MidiTrack::omniModeOff(unsigned int delta) {
    addData(uintToVLQ(delta));
    addData({CHANNEL_MODE, 0x7C, 0x00});
}

//receiver should respong to voice messages on any MIDI channel
void MidiTrack::omniModeOn(unsigned int delta) {
    addData(uintToVLQ(delta));
    addData({CHANNEL_MODE, 0x7D, 0x00});
}

//puts the receiver into monophonic mode
void MidiTrack::monoModeOn(unsigned int delta, byte channels) {
    addData(uintToVLQ(delta));
    addData({CHANNEL_MODE, 0x7E, channels});
}

//puts the receiver into polyphonic mode
void MidiTrack::polyModeOn(unsigned int delta) {
    addData(uintToVLQ(delta));
    addData({CHANNEL_MODE, 0x7F, 0x00});
}


//--MIDI FILE IMPLEMENTATION


//writes an arbitrary chunk to the midi file
void MidiFile::writeChunk(MidiChunk chunk) {

    /* get info from chunk only once
     * which should lead to less function calls
     * which in turn may lead to better time efficiency
     * at the expense of a little extra memory (I think)
     */
    std::vector<byte> typeAndLength = chunk.getType();
    uint32_t length = chunk.getLength();
    std::vector<byte> data = chunk.getData();

    //combine type and length so I need only iterate through two vectors (neater code)
    typeAndLength.push_back(getNthByte32(length,3));
    typeAndLength.push_back(getNthByte32(length,2));
    typeAndLength.push_back(getNthByte32(length,1));
    typeAndLength.push_back(getNthByte32(length,0));

    //write vectors to my file
    std::ostream_iterator<byte> writeBytes(midiFile);
    std::copy(std::begin(typeAndLength), std::end(typeAndLength), writeBytes);
    std::copy(std::begin(data), std::end(data), writeBytes);
}

//needs a name for the file, and a header file and an initial track (can have more)
MidiFile::MidiFile(std::string fileName, MidiHeader hd, MidiTrack trk) {
    midiFile.open(fileName, std::ofstream::out | std::ofstream::trunc);
    writeChunk(std::move(hd));
    writeChunk(std::move(trk));
}

/* writes a track to the midi file,
 * just a call to writeChunk, but enforced to be a track
 */
void MidiFile::addTrack(MidiTrack trk) {
    writeChunk(std::move(trk));
}

void MidiFile::closeFile() {
    midiFile.close();
}

//--IMPLEMENTATION OF UTILITY FUNCTIONS--

VLQ uintToVLQ(unsigned int num) {
    VLQ converted;
    if(num == 0) {
        converted.push_back(0x00);
    } else {
        converted.push_back((byte)(num & 0x7F)); //this is the end byte, so has bit 7 set to 0
        num = num >> 7; //move along 7 bits
        while(num != 0) {
            converted.push_back((byte)((num & 0x7F)|0x80)); // add 1 as msb in byte as not at the end of VLQ value
            num = num >> 7; //shift right another 7 bits
        }
    }
    std::reverse(std::begin(converted),std::end(converted)); //reversing the vector as i started with lsb first
    return converted;

}

/* credit for this function goes to Nawaz at:
 * https://stackoverflow.com/questions/7378087/how-to-efficiently-copy-a-stdstring-into-a-vector
 */
std::vector<byte> fromString(std::string str) {
    std::vector<byte> byteForm (std::begin(str),std::end(str));
    byteForm.push_back('\0'); //since being written as raw bytes, adding a \0 is probably a sensible ides
    return byteForm;
}

/*credit for how to get byte n goes to oon at:
 * https://stackoverflow.com/questions/7277926/bitwise-operator-to-get-byte-from-32-bits
 */
byte getNthByte16(uint16_t num, unsigned int n) {
    if(n <= 1) {
        return (byte)((num >> (n << 3)) & 0xFF); //no real need for cast, but suppresses warning
    } else {
        throw "Invalid value of n for 16 bit value (0-1 valid)";
    }
}

byte getNthByte32(uint32_t num, unsigned int n) {
    if(n <= 3) {
        return (byte)((num >> (n << 3)) & 0xFF); //cast used for removing warnings
    } else {
        throw "Invalid value of n for 32 bit value (0-3 valid)";
    }
}
