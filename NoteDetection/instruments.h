#include <vector>

class Instrument 
{
public:
    virtual float getMinFreq();
    virtual float getMaxFreq();
    virtual const int* getHarmonics();
};

class Bass : public Instrument 
{
public:
    virtual float getMinFreq() override { return 21.1f; };
    virtual float getMaxFreq() override { return 5200.0f; };
};

class Guitar : public Instrument 
{
public:
    virtual float getMinFreq() override { return 82.41f; };
    virtual float getMaxFreq() override { return 5200.0f; };
};

class Piano : public Instrument 
{
public:
    virtual float getMinFreq() override { return 15.0f; };
    virtual float getMaxFreq() override { return 16000.0f; };
};

class Violin : public Instrument 
{
public:
    virtual float getMinFreq() override { return 196.0f; };
    virtual float getMaxFreq() override { return 2637.0f; };
};

class Voice : public Instrument 
{
public:
    virtual float getMinFreq() override { return 21.1f; };
    virtual float getMaxFreq() override { return 880.0; /*A5*/ };
};