/*
  ==============================================================================

    This file contains the basic framework code for a JUCE plugin editor.

  ==============================================================================
*/

#pragma once

#include <JuceHeader.h>
#include "PluginProcessor.h"

//==============================================================================
/**
*/
class SimpleTunerAudioProcessorEditor  : public juce::AudioProcessorEditor
{
public:
    SimpleTunerAudioProcessorEditor (SimpleTunerAudioProcessor&);
    ~SimpleTunerAudioProcessorEditor() override;

    //==============================================================================
    void paint (juce::Graphics&) override;
    void resized() override;

private:
    // This reference is provided as a quick way for your editor to
    // access the processor object that created it.
    SimpleTunerAudioProcessor& audioProcessor;
    
    // adapted from tutorial https://docs.juce.com/master/tutorial_label.html
    juce::Label titleLabel;
    juce::Label inputLabel;
    juce::Label inputText;
    juce::Label uppercaseLabel;
    juce::Label uppercaseText;
    juce::Label trailLabel;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (SimpleTunerAudioProcessorEditor)
};
