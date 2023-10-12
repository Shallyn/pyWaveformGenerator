# pyWaveformGenerator
calculate GW waveform

Installation && compile
```
./mkconf.sh
./configure
make
```
sometimes the lib cannot use in Apple core macbook.


Run in python
```
from pyWaveformGenerator.waveform import calculate_waveform
waveform, dynamics = calculate_waveform((10, 10, 
	0, 0, 0, 
	0, 0, 0, 
	0, 100, 
	0, 0, 0, 0), 20)
```
see more details and description in waveform.py
