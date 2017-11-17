import pyaudio
import wave
import audioop
from spinn_front_end_common.utilities.connections import LiveEventConnection



FORMAT = pyaudio.paFloat32
CHANNELS = 1
RATE = 22050
CHUNK = 1024
RECORD_SECONDS = 1
WAVE_OUTPUT_FILENAME = "file.wav"

audio = pyaudio.PyAudio()

#Live event connection
live_samples = LiveEventConnection()

# start Recording
stream = audio.open(format=FORMAT, channels=CHANNELS,
                    rate=RATE, input=True,
                    frames_per_buffer=CHUNK)
print "recording..."
frames = []
maximum=0

for i in range(0, int(RATE / CHUNK * RECORD_SECONDS)):
    data = stream.read(CHUNK)

    #convert to SPL

    #send the data chunk to the live event connection
    #TODO: check if MC payload can be sent from live event connection


    mx=audioop.max(data, 2)
    if mx>maximum:
        maximum=mx
    frames.append(data)
print "finished recording"
print maximum
# stop Recording
stream.stop_stream()
stream.close()
audio.terminate()

waveFile = wave.open(WAVE_OUTPUT_FILENAME, 'wb')
waveFile.setnchannels(CHANNELS)
waveFile.setsampwidth(audio.get_sample_size(FORMAT))
waveFile.setframerate(RATE)
waveFile.writeframes(b''.join(frames))
waveFile.close()