/*
  ==============================================================================

    This file contains the basic framework code for a JUCE plugin editor.

  ==============================================================================
*/

#include "PluginProcessor.h"
#include "PluginEditor.h"

//==============================================================================
SimpleTunerAudioProcessorEditor::SimpleTunerAudioProcessorEditor (SimpleTunerAudioProcessor& p)
    : AudioProcessorEditor (&p), audioProcessor (p)
{
    // Make sure that before the constructor has finished, you've set the
    // editor's size to whatever you need it to be.
    // component?
   
    // loosely copied a tutoria for gui components trying to figure out how it works : https://docs.juce.com/master/tutorial_label.html

    addAndMakeVisible(titleLabel);
    titleLabel.setFont(juce::Font(16.0f, juce::Font::bold));
    titleLabel.setText("Simple Chromatic Tuner", juce::dontSendNotification);
    titleLabel.setColour(juce::Label::textColourId, juce::Colours::lightgreen);
    titleLabel.setJustificationType(juce::Justification::centred);

    addAndMakeVisible(inputLabel);
    inputLabel.setText("A4 Reference Pitch:", juce::dontSendNotification);
    inputLabel.attachToComponent(&inputText, true);
    inputLabel.setColour(juce::Label::textColourId, juce::Colours::orange);
    inputLabel.setJustificationType(juce::Justification::right);

    addAndMakeVisible(uppercaseLabel);
    uppercaseLabel.setText("Uppercase:", juce::dontSendNotification);
    uppercaseLabel.attachToComponent(&uppercaseText, true);
    uppercaseLabel.setColour(juce::Label::textColourId, juce::Colours::orange);
    uppercaseLabel.setJustificationType(juce::Justification::right);

    addAndMakeVisible(inputText);
    inputText.setEditable(true);
    inputText.setColour(juce::Label::backgroundColourId, juce::Colours::darkblue);
    inputText.setText("440",juce::dontSendNotification);
    inputText.onTextChange = [this] { uppercaseText.setText(inputText.getText().toUpperCase(), juce::dontSendNotification); };
    inputText.attachToComponent(&trailLabel, true);

    addAndMakeVisible(uppercaseText);
    uppercaseText.setColour(juce::Label::backgroundColourId, juce::Colours::darkblue);
    uppercaseText.setText(inputText.getText().toUpperCase(), juce::dontSendNotification);

    addAndMakeVisible(trailLabel);
    trailLabel.setText("Hz", juce::dontSendNotification);
    //inputLabel.attachToComponent(&inputText, true);
    trailLabel.setColour(juce::Label::textColourId, juce::Colours::orange);
    trailLabel.setJustificationType(juce::Justification::left);

    setResizable(true, true);
    setSize(320, 200);

}

SimpleTunerAudioProcessorEditor::~SimpleTunerAudioProcessorEditor()
{
}

//==============================================================================
void SimpleTunerAudioProcessorEditor::paint (juce::Graphics& g)
{
    // (Our component is opaque, so we must completely fill the background with a solid colour)
    g.fillAll (getLookAndFeel().findColour (juce::ResizableWindow::backgroundColourId));

    g.setColour(juce::Colours::goldenrod);
    g.setFont (15.0f);
    g.drawFittedText ("W = " + std::to_string(getWidth()) + " H = " + std::to_string(getHeight()), getLocalBounds(), juce::Justification::centredBottom, 1);

     

}

void SimpleTunerAudioProcessorEditor::resized()
{
    // This is generally where you'll want to lay out the positions of any
    // subcomponents in your editor..

    titleLabel.setBounds(10, 10, getWidth() - 20, 30);
    trailLabel.setBounds(getWidth() - 30, 50, getWidth(), 20); // lame (need a better way of layouting)
    uppercaseText.setBounds(100, 80, getWidth() - 110, 20);
}
