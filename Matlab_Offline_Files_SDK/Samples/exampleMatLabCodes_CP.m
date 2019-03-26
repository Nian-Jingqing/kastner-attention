MATLAB is selecting SOFTWARE OPENGL rendering.

                            < M A T L A B (R) >
                  Copyright 1984-2017 The MathWorks, Inc.
                   R2017b (9.3.0.713579) 64-bit (glnxa64)
                             September 14, 2017

 
To get started, type one of these: helpwin, helpdesk, or demo.
For product information, visit www.mathworks.com.
 
>> addpath('/snel/home/fzhu23/emacs/matlab-emacs-src/toolbox','-begin'); rehash; emacsinit('emacsclient -n');
>> 
>> ls /mnt/scratch/feng
RP_06182018_SPKC.pl2

>> pwd

ans =

    '/snel/home/fzhu23'

>> cd bin
>> 
>> cd analysis_tools/
>> ls
+Continuous	    +GLM		Matlab Offline Files SDK  README.md
+DataPreprocessing  +GPFA		+Movement		  +Run
+Datasets	    gpfa_utils		+Plot			  +Tensor
+Examples	    +jPCA_with_CP_mods	+R			  +Utils

>> 
>> cd Matlab Offline Files SDK/
Error using cd
Too many input arguments.
 
>> cd('Matlab Offline Files SDK/')
>> ls
Change Log.txt			     PL2ReadNextDataBlock.m
ddt.m				     PL2StartStopTs.m
ddt_v.m				     PL2TsBySource.m
ddt_write_v.m			     PL2Ts.m
internalPL2Ad.p			     PL2WavesBySource.m
internalPL2AdSpan.p		     PL2Waves.m
internalPL2AdTimeSpan.p		     plx_adchan_freqs.m
internalPL2EventTs.p		     plx_adchan_gains.m
internalPL2GetChannel.p		     plx_ad_chanmap.m
internalPL2MergeFragments.p	     plx_adchan_names.m
internalPL2ReadFileIndex.p	     plx_adchan_samplecounts.m
internalPL2ReadFirstDataBlock.p      plx_ad_gap_info.m
internalPL2ReadNextDataBlock.p	     plx_ad_info.m
internalPL2ReadPdp.p		     plx_ad.m
internalPL2ResolveChannelBySource.p  plx_ad_resolve_channel.m
internalPL2ResolveChannel.p	     plx_ad_span.m
internalPL2ResolveFilename.p	     plx_ad_span_v.m
internalPL2ResolveFilenamePlx.p      plx_ad_v.m
internalPL2StartStopTs.p	     plx_chan_filters.m
internalPL2TrimNames.p		     plx_chan_gains.m
internalPL2Ts.p			     plx_chanmap.m
internalPL2VerifyRecord.p	     plx_chan_names.m
internalPL2Waves.p		     plx_chan_thresholds.m
Matlab Offline Files SDK.pdf	     plx_close.m
mexPlex				     plx_event_chanmap.m
mexPlex.mexglx			     plx_event_names.m
mexPlex.mexw32			     plx_event_resolve_channel.m
mexPlex.mexw64			     plx_event_ts.m
PL2AdBySource.m			     plx_info.m
PL2Ad.m				     plx_information.m
PL2AdSpanBySource.m		     plx_mexplex_version.m
PL2AdSpan.m			     plx_resolve_channel.m
PL2AdTimeSpanBySource.m		     plx_spike_info.m
PL2AdTimeSpan.m			     plx_ts.m
PL2EventTsBySource.m		     plx_vt_interpret.m
PL2EventTs.m			     plx_waves.m
PL2GetFileIndex.m		     plx_waves_v.m
PL2Print.m			     readall.m
PL2ReadFileIndex.m		     Samples
PL2ReadFirstDataBlock.m

>> 
>> cd Samples/
>> readall
Opened File Name: /mnt/scratch/feng/RP_06182018_SPKC.pl2
Version: 3
Frequency : 40000
Comment : 
Date/Time :  6/18/2018 13:37: 1
Duration : 11268.9759
Num Pts Per Wave : 56
Num Pts Pre-Threshold : 16
loading channel 129 / 384
Elapsed time is 105.119957 seconds.
loading channel 130 / 384
Elapsed time is 36.692190 seconds.
loading channel 131 / 384
Elapsed time is 75.119971 seconds.
loading channel 132 / 384
  C-c C-cOperation terminated by user during internalPL2Ad


In PL2Ad (line 32)
[ad] = internalPL2Ad(filename, channel);

In plx_ad (line 43)
    pl2ad = PL2Ad(filename, pl2Channel);

In readall (line 89)
            [adfreq, nad, tsad, fnad, allad{ich+1}] = plx_ad(OpenedFileName,
            ich);
 
>> whos
  Name                    Size                   Bytes  Class     Attributes

  Comment                 0x0                        0  char                
  DateTime                1x19                      38  char                
  Duration                1x1                        8  double              
  Freq                    1x1                        8  double              
  NPW                     1x1                        8  double              
  OpenedFileName          1x38                      76  char                
  PreThresh               1x1                        8  double              
  SlowADResBits           1x1                        8  double              
  SlowPeakV               1x1                        8  double              
  SpikeADResBits          1x1                        8  double              
  SpikePeakV              1x1                        8  double              
  StartingFileName        1x38                      76  char                
  Trodalness              1x1                        8  double              
  Version                 1x1                        8  double              
  adfreq                  1x1                        8  double              
  allad                   1x384            10818211824  cell                
  allts                  27x129                  27864  cell                
  ans                     1x17                      34  char                
  evcounts                1x43                     344  double              
  fnad                    1x1                        8  double              
  ich                     1x1                        8  double              
  iunit                   1x1                        8  double              
  nad                     1x1                        8  double              
  nchannels1              1x1                        8  double              
  nslowchannels           1x1                        8  double              
  nspk                    1x1                        8  double              
  numads                  1x1                        8  double              
  nunits1                 1x1                        8  double              
  slowcounts              1x384                   3072  double              
  spk_filters           128x1                     1024  double              
  spk_gains             128x1                     1024  double              
  spk_names             128x6                     1536  char                
  spk_threshs           128x1                     1024  double              
  tsad                    1x1                        8  double              
  tscounts               27x129                  27864  double              
  u                       1x1                        8  double              
  wfcounts               27x129                  27864  double              

>> size(allad)

ans =

     1   384

>> allad

allad =

  1x384 cell array

  Columns 1 through 4

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 5 through 8

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 9 through 12

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 13 through 16

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 17 through 20

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 21 through 24

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 25 through 28

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 29 through 32

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 33 through 36

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 37 through 40

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 41 through 44

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 45 through 48

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 49 through 52

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 53 through 56

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 57 through 60

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 61 through 64

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 65 through 68

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 69 through 72

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 73 through 76

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 77 through 80

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 81 through 84

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 85 through 88

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 89 through 92

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 93 through 96

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 97 through 100

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 101 through 104

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 105 through 108

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 109 through 112

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 113 through 116

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 117 through 120

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 121 through 124

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 125 through 128

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 129 through 131

    {450758685x1 double}    {450758685x1 double}    {450758685x1 double}

  Columns 132 through 135

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 136 through 139

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 140 through 143

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 144 through 147

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 148 through 151

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 152 through 155

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 156 through 159

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 160 through 163

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 164 through 167

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 168 through 171

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 172 through 175

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 176 through 179

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 180 through 183

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 184 through 187

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 188 through 191

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 192 through 195

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 196 through 199

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 200 through 203

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 204 through 207

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 208 through 211

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 212 through 215

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 216 through 219

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 220 through 223

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 224 through 227

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 228 through 231

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 232 through 235

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 236 through 239

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 240 through 243

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 244 through 247

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 248 through 251

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 252 through 255

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 256 through 259

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 260 through 263

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 264 through 267

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 268 through 271

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 272 through 275

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 276 through 279

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 280 through 283

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 284 through 287

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 288 through 291

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 292 through 295

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 296 through 299

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 300 through 303

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 304 through 307

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 308 through 311

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 312 through 315

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 316 through 319

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 320 through 323

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 324 through 327

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 328 through 331

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 332 through 335

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 336 through 339

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 340 through 343

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 344 through 347

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 348 through 351

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 352 through 355

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 356 through 359

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 360 through 363

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 364 through 367

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 368 through 371

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 372 through 375

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 376 through 379

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 380 through 383

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Column 384

    {0x0 double}

>> tmp = horzcat(allad {[129 130 131]} );
size(tmp)
>>   C-c C-c  C-c C-c
>> 
>> 
>> size(tmp)

ans =

   450758685           3

>> plot(tmp(1:10000,:))
Warning: MATLAB has disabled some advanced graphics rendering features by
switching to software OpenGL. For more information, click here. 
>> tmp2 = tmp(1:100000,:) - mean(tmp(1:100000,:),2);
>> plot(tmp2)
>> help plx_ad
  plx_ad(filename, channel): read all a/d data for the specified channel from a .plx or .pl2 file
               returns raw a/d values
 
  [adfreq, n, ts, fn, ad] = plx_ad(filename, channel)
  [adfreq, n, ts, fn, ad] = plx_ad(filename, 0)
  [adfreq, n, ts, fn, ad] = plx_ad(filename, 'FP01')
 
  INPUT:
    filename - if empty string, will use File Open dialog
    channel - 0-based channel number or channel name
 
            a/d data come in fragments. Each fragment has a timestamp
            and a number of a/d data points. The timestamp corresponds to
            the time of recording of the first a/d value in this fragment.
            All the data values stored in the vector ad.
  
  OUTPUT:
    adfreq - digitization frequency for this channel
    n - total number of data points 
    ts - array of fragment timestamps (one timestamp per fragment, in seconds)
    fn - number of data points in each fragment
    ad - array of raw a/d values

>> readall
Opened File Name: /mnt/scratch/feng/RP_06182018_SPKC.pl2
Version: 3
Frequency : 40000
Comment : 
Date/Time :  6/18/2018 13:37: 1
Duration : 11268.9759
Num Pts Per Wave : 56
Num Pts Pre-Threshold : 16
loading channel 129 / 384
Error using plx_ad_span
Too many output arguments.

Error in readall (line 90)
            [adfreq, nad, tsad, fnad, allad{ich+1}] = ...
 
>> readall
Opened File Name: /mnt/scratch/feng/RP_06182018_SPKC.pl2
Version: 3
Frequency : 40000
Comment : 
Date/Time :  6/18/2018 13:37: 1
Duration : 11268.9759
Num Pts Per Wave : 56
Num Pts Pre-Threshold : 16
loading channel 129 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.024799 seconds.
loading channel 130 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.001341 seconds.
loading channel 131 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.001201 seconds.
loading channel 132 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.001559 seconds.
loading channel 133 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.005778 seconds.
loading channel 134 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000693 seconds.
loading channel 135 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000560 seconds.
loading channel 136 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000544 seconds.
loading channel 137 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000748 seconds.
loading channel 138 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000532 seconds.
loading channel 139 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000527 seconds.
loading channel 140 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000513 seconds.
loading channel 141 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000523 seconds.
loading channel 142 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000522 seconds.
loading channel 143 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000509 seconds.
loading channel 144 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000517 seconds.
loading channel 145 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000506 seconds.
loading channel 146 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000518 seconds.
loading channel 147 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000513 seconds.
loading channel 148 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000507 seconds.
loading channel 149 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000518 seconds.
loading channel 150 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000507 seconds.
loading channel 151 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000516 seconds.
loading channel 152 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000518 seconds.
loading channel 153 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000511 seconds.
loading channel 154 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000516 seconds.
loading channel 155 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000507 seconds.
loading channel 156 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000513 seconds.
loading channel 157 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000519 seconds.
loading channel 158 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000525 seconds.
loading channel 159 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000532 seconds.
loading channel 160 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000509 seconds.
loading channel 161 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000518 seconds.
loading channel 162 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000517 seconds.
loading channel 163 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000510 seconds.
loading channel 164 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000519 seconds.
loading channel 165 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000508 seconds.
loading channel 166 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000516 seconds.
loading channel 167 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000526 seconds.
loading channel 168 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000510 seconds.
loading channel 169 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000515 seconds.
loading channel 170 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000506 seconds.
loading channel 171 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000513 seconds.
loading channel 172 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000516 seconds.
loading channel 173 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000507 seconds.
loading channel 174 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000518 seconds.
loading channel 175 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000515 seconds.
loading channel 176 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000517 seconds.
loading channel 177 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000512 seconds.
loading channel 178 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000505 seconds.
loading channel 179 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000515 seconds.
loading channel 180 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000510 seconds.
loading channel 181 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000515 seconds.
loading channel 182 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000515 seconds.
loading channel 183 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000510 seconds.
loading channel 184 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000518 seconds.
loading channel 185 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000506 seconds.
loading channel 186 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000515 seconds.
loading channel 187 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000515 seconds.
loading channel 188 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000508 seconds.
loading channel 189 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000518 seconds.
loading channel 190 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000504 seconds.
loading channel 191 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000514 seconds.
loading channel 192 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000515 seconds.
loading channel 193 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000507 seconds.
loading channel 194 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000519 seconds.
loading channel 195 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000506 seconds.
loading channel 196 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000515 seconds.
loading channel 197 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000512 seconds.
loading channel 198 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000506 seconds.
loading channel 199 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000514 seconds.
loading channel 200 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000505 seconds.
loading channel 201 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000515 seconds.
loading channel 202 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000505 seconds.
loading channel 203 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000504 seconds.
loading channel 204 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000515 seconds.
loading channel 205 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000515 seconds.
loading channel 206 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000517 seconds.
loading channel 207 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000507 seconds.
loading channel 208 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000505 seconds.
loading channel 209 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000512 seconds.
loading channel 210 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000502 seconds.
loading channel 211 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000513 seconds.
loading channel 212 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000502 seconds.
loading channel 213 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000505 seconds.
loading channel 214 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000513 seconds.
loading channel 215 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000504 seconds.
loading channel 216 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000518 seconds.
loading channel 217 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000507 seconds.
loading channel 218 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000510 seconds.
loading channel 219 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000510 seconds.
loading channel 220 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000505 seconds.
loading channel 221 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000516 seconds.
loading channel 222 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000503 seconds.
loading channel 223 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000513 seconds.
loading channel 224 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000514 seconds.
loading channel 225 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000505 seconds.
loading channel 226 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000513 seconds.
loading channel 227 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000505 seconds.
loading channel 228 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000513 seconds.
loading channel 229 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000507 seconds.
loading channel 230 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000506 seconds.
loading channel 231 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000513 seconds.
loading channel 232 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000500 seconds.
loading channel 233 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000513 seconds.
loading channel 234 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000502 seconds.
loading channel 235 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000503 seconds.
loading channel 236 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000512 seconds.
loading channel 237 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000504 seconds.
loading channel 238 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000512 seconds.
loading channel 239 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000503 seconds.
loading channel 240 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000504 seconds.
loading channel 241 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000512 seconds.
loading channel 242 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000504 seconds.
loading channel 243 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000514 seconds.
loading channel 244 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000501 seconds.
loading channel 245 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000511 seconds.
loading channel 246 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000510 seconds.
loading channel 247 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000504 seconds.
loading channel 248 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000519 seconds.
loading channel 249 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000503 seconds.
loading channel 250 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000509 seconds.
loading channel 251 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000509 seconds.
loading channel 252 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000501 seconds.
loading channel 253 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000513 seconds.
loading channel 254 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000501 seconds.
loading channel 255 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000509 seconds.
loading channel 256 / 384

 PL2AdSpan: startCount should be positive.
Elapsed time is 0.000501 seconds.
>> allad

allad =

  1x384 cell array

  Columns 1 through 4

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 5 through 8

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 9 through 12

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 13 through 16

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 17 through 20

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 21 through 24

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 25 through 28

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 29 through 32

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 33 through 36

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 37 through 40

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 41 through 44

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 45 through 48

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 49 through 52

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 53 through 56

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 57 through 60

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 61 through 64

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 65 through 68

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 69 through 72

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 73 through 76

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 77 through 80

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 81 through 84

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 85 through 88

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 89 through 92

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 93 through 96

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 97 through 100

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 101 through 104

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 105 through 108

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 109 through 112

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 113 through 116

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 117 through 120

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 121 through 124

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 125 through 129

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}    {[-1]}

  Columns 130 through 136

    {[-1]}    {[-1]}    {[-1]}    {[-1]}    {[-1]}    {[-1]}    {[-1]}

  Columns 137 through 143

    {[-1]}    {[-1]}    {[-1]}    {[-1]}    {[-1]}    {[-1]}    {[-1]}

  Columns 144 through 150

    {[-1]}    {[-1]}    {[-1]}    {[-1]}    {[-1]}    {[-1]}    {[-1]}

  Columns 151 through 157

    {[-1]}    {[-1]}    {[-1]}    {[-1]}    {[-1]}    {[-1]}    {[-1]}

  Columns 158 through 164

    {[-1]}    {[-1]}    {[-1]}    {[-1]}    {[-1]}    {[-1]}    {[-1]}

  Columns 165 through 171

    {[-1]}    {[-1]}    {[-1]}    {[-1]}    {[-1]}    {[-1]}    {[-1]}

  Columns 172 through 178

    {[-1]}    {[-1]}    {[-1]}    {[-1]}    {[-1]}    {[-1]}    {[-1]}

  Columns 179 through 185

    {[-1]}    {[-1]}    {[-1]}    {[-1]}    {[-1]}    {[-1]}    {[-1]}

  Columns 186 through 192

    {[-1]}    {[-1]}    {[-1]}    {[-1]}    {[-1]}    {[-1]}    {[-1]}

  Columns 193 through 199

    {[-1]}    {[-1]}    {[-1]}    {[-1]}    {[-1]}    {[-1]}    {[-1]}

  Columns 200 through 206

    {[-1]}    {[-1]}    {[-1]}    {[-1]}    {[-1]}    {[-1]}    {[-1]}

  Columns 207 through 213

    {[-1]}    {[-1]}    {[-1]}    {[-1]}    {[-1]}    {[-1]}    {[-1]}

  Columns 214 through 220

    {[-1]}    {[-1]}    {[-1]}    {[-1]}    {[-1]}    {[-1]}    {[-1]}

  Columns 221 through 227

    {[-1]}    {[-1]}    {[-1]}    {[-1]}    {[-1]}    {[-1]}    {[-1]}

  Columns 228 through 234

    {[-1]}    {[-1]}    {[-1]}    {[-1]}    {[-1]}    {[-1]}    {[-1]}

  Columns 235 through 241

    {[-1]}    {[-1]}    {[-1]}    {[-1]}    {[-1]}    {[-1]}    {[-1]}

  Columns 242 through 248

    {[-1]}    {[-1]}    {[-1]}    {[-1]}    {[-1]}    {[-1]}    {[-1]}

  Columns 249 through 255

    {[-1]}    {[-1]}    {[-1]}    {[-1]}    {[-1]}    {[-1]}    {[-1]}

  Columns 256 through 260

    {[-1]}    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 261 through 264

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 265 through 268


re
  Columns 269 through 272

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 273 through 276

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 277 through 280

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 281 through 284

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 285 through 288

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 289 through 292

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 293 through 296

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 297 through 300

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 301 through 304

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 305 through 308

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 309 through 312

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 313 through 316

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 317 through 320

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 321 through 324

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 325 through 328

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 329 through 332

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 333 through 336

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 337 through 340

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 341 through 344

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 345 through 348

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 349 through 352

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 353 through 356

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 357 through 360

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 361 through 364

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 365 through 368

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 369 through 372

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 373 through 376

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 377 through 380

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 381 through 384

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

>> readall
Opened File Name: /mnt/scratch/feng/RP_06182018_SPKC.pl2
Version: 3
Frequency : 40000
Comment : 
Date/Time :  6/18/2018 13:37: 1
Duration : 11268.9759
Num Pts Per Wave : 56
Num Pts Pre-Threshold : 16
loading channel 129 / 384
Elapsed time is 0.198141 seconds.
loading channel 130 / 384
Elapsed time is 0.032430 seconds.
loading channel 131 / 384
Elapsed time is 0.026763 seconds.
loading channel 132 / 384
Elapsed time is 0.029479 seconds.
loading channel 133 / 384
Elapsed time is 0.037125 seconds.
loading channel 134 / 384
Elapsed time is 0.031056 seconds.
loading channel 135 / 384
Elapsed time is 0.041658 seconds.
loading channel 136 / 384
Elapsed time is 0.027617 seconds.
loading channel 137 / 384
Elapsed time is 0.027795 seconds.
loading channel 138 / 384
Elapsed time is 0.028321 seconds.
loading channel 139 / 384
Elapsed time is 0.028610 seconds.
loading channel 140 / 384
Elapsed time is 0.030638 seconds.
loading channel 141 / 384
Elapsed time is 0.029391 seconds.
loading channel 142 / 384
Elapsed time is 0.028181 seconds.
loading channel 143 / 384
Elapsed time is 0.027487 seconds.
loading channel 144 / 384
Elapsed time is 0.028037 seconds.
loading channel 145 / 384
Elapsed time is 0.030529 seconds.
loading channel 146 / 384
Elapsed time is 0.028841 seconds.
loading channel 147 / 384
Elapsed time is 0.030611 seconds.
loading channel 148 / 384
Elapsed time is 0.031844 seconds.
loading channel 149 / 384
Elapsed time is 0.031404 seconds.
loading channel 150 / 384
Elapsed time is 0.028302 seconds.
loading channel 151 / 384
Elapsed time is 0.027182 seconds.
loading channel 152 / 384
Elapsed time is 0.026859 seconds.
loading channel 153 / 384
Elapsed time is 0.027530 seconds.
loading channel 154 / 384
Elapsed time is 0.027548 seconds.
loading channel 155 / 384
Elapsed time is 0.034789 seconds.
loading channel 156 / 384
Elapsed time is 0.031226 seconds.
loading channel 157 / 384
Elapsed time is 0.029971 seconds.
loading channel 158 / 384
Elapsed time is 0.028031 seconds.
loading channel 159 / 384
Elapsed time is 0.027183 seconds.
loading channel 160 / 384
Elapsed time is 0.029824 seconds.
loading channel 161 / 384
Elapsed time is 0.031837 seconds.
loading channel 162 / 384
Elapsed time is 0.028451 seconds.
loading channel 163 / 384
Elapsed time is 0.027929 seconds.
loading channel 164 / 384
Elapsed time is 0.028774 seconds.
loading channel 165 / 384
Elapsed time is 0.030713 seconds.
loading channel 166 / 384
Elapsed time is 0.030988 seconds.
loading channel 167 / 384
Elapsed time is 0.027905 seconds.
loading channel 168 / 384
Elapsed time is 0.027015 seconds.
loading channel 169 / 384
Elapsed time is 0.027067 seconds.
loading channel 170 / 384
Elapsed time is 0.027419 seconds.
loading channel 171 / 384
Elapsed time is 0.027590 seconds.
loading channel 172 / 384
Elapsed time is 0.028816 seconds.
loading channel 173 / 384
Elapsed time is 0.029036 seconds.
loading channel 174 / 384
Elapsed time is 0.028692 seconds.
loading channel 175 / 384
Elapsed time is 0.027838 seconds.
loading channel 176 / 384
Elapsed time is 0.027184 seconds.
loading channel 177 / 384
Elapsed time is 0.027718 seconds.
loading channel 178 / 384
Elapsed time is 0.029617 seconds.
loading channel 179 / 384
Elapsed time is 0.032775 seconds.
loading channel 180 / 384
Elapsed time is 0.032104 seconds.
loading channel 181 / 384
Elapsed time is 0.031292 seconds.
loading channel 182 / 384
Elapsed time is 0.040134 seconds.
loading channel 183 / 384
Elapsed time is 0.030279 seconds.
loading channel 184 / 384
Elapsed time is 0.028503 seconds.
loading channel 185 / 384
Elapsed time is 0.026015 seconds.
loading channel 186 / 384
Elapsed time is 0.026111 seconds.
loading channel 187 / 384
Elapsed time is 0.026899 seconds.
loading channel 188 / 384
Elapsed time is 0.031903 seconds.
loading channel 189 / 384
Elapsed time is 0.030877 seconds.
loading channel 190 / 384
Elapsed time is 0.027389 seconds.
loading channel 191 / 384
Elapsed time is 0.027538 seconds.
loading channel 192 / 384
Elapsed time is 0.028667 seconds.
loading channel 193 / 384
Elapsed time is 0.028420 seconds.
loading channel 194 / 384
Elapsed time is 0.029309 seconds.
loading channel 195 / 384
Elapsed time is 0.029087 seconds.
loading channel 196 / 384
Elapsed time is 0.032196 seconds.
loading channel 197 / 384
Elapsed time is 0.032249 seconds.
loading channel 198 / 384
Elapsed time is 0.028085 seconds.
loading channel 199 / 384
Elapsed time is 0.037699 seconds.
loading channel 200 / 384
Elapsed time is 0.027952 seconds.
loading channel 201 / 384
Elapsed time is 0.026309 seconds.
loading channel 202 / 384
Elapsed time is 0.026196 seconds.
loading channel 203 / 384
Elapsed time is 0.026443 seconds.
loading channel 204 / 384
Elapsed time is 0.027322 seconds.
loading channel 205 / 384
Elapsed time is 0.027460 seconds.
loading channel 206 / 384
Elapsed time is 0.026436 seconds.
loading channel 207 / 384
Elapsed time is 0.025737 seconds.
loading channel 208 / 384
Elapsed time is 0.025980 seconds.
loading channel 209 / 384
Elapsed time is 0.026828 seconds.
loading channel 210 / 384
Elapsed time is 0.026120 seconds.
loading channel 211 / 384
Elapsed time is 0.027711 seconds.
loading channel 212 / 384
Elapsed time is 0.028952 seconds.
loading channel 213 / 384
Elapsed time is 0.029219 seconds.
loading channel 214 / 384
Elapsed time is 0.027253 seconds.
loading channel 215 / 384
Elapsed time is 0.028934 seconds.
loading channel 216 / 384
Elapsed time is 0.028455 seconds.
loading channel 217 / 384
Elapsed time is 0.035074 seconds.
loading channel 218 / 384
Elapsed time is 0.029483 seconds.
loading channel 219 / 384
Elapsed time is 0.028979 seconds.
loading channel 220 / 384
Elapsed time is 0.033508 seconds.
loading channel 221 / 384
Elapsed time is 0.033163 seconds.
loading channel 222 / 384
Elapsed time is 0.029374 seconds.
loading channel 223 / 384
Elapsed time is 0.027675 seconds.
loading channel 224 / 384
Elapsed time is 0.027786 seconds.
loading channel 225 / 384
Elapsed time is 0.028374 seconds.
loading channel 226 / 384
Elapsed time is 0.029123 seconds.
loading channel 227 / 384
Elapsed time is 0.032401 seconds.
loading channel 228 / 384
Elapsed time is 0.034582 seconds.
loading channel 229 / 384
Elapsed time is 0.029119 seconds.
loading channel 230 / 384
Elapsed time is 0.027842 seconds.
loading channel 231 / 384
Elapsed time is 0.029839 seconds.
loading channel 232 / 384
Elapsed time is 0.028556 seconds.
loading channel 233 / 384
Elapsed time is 0.027443 seconds.
loading channel 234 / 384
Elapsed time is 0.027970 seconds.
loading channel 235 / 384
Elapsed time is 0.027556 seconds.
loading channel 236 / 384
Elapsed time is 0.028793 seconds.
loading channel 237 / 384
Elapsed time is 0.028856 seconds.
loading channel 238 / 384
Elapsed time is 0.028075 seconds.
loading channel 239 / 384
Elapsed time is 0.028104 seconds.
loading channel 240 / 384
Elapsed time is 0.029452 seconds.
loading channel 241 / 384
Elapsed time is 0.031601 seconds.
loading channel 242 / 384
Elapsed time is 0.033972 seconds.
loading channel 243 / 384
Elapsed time is 0.030409 seconds.
loading channel 244 / 384
Elapsed time is 0.032072 seconds.
loading channel 245 / 384
Elapsed time is 0.033503 seconds.
loading channel 246 / 384
Elapsed time is 0.029379 seconds.
loading channel 247 / 384
Elapsed time is 0.030961 seconds.
loading channel 248 / 384
Elapsed time is 0.029158 seconds.
loading channel 249 / 384
Elapsed time is 0.026178 seconds.
loading channel 250 / 384
Elapsed time is 0.026499 seconds.
loading channel 251 / 384
Elapsed time is 0.037027 seconds.
loading channel 252 / 384
Elapsed time is 0.033237 seconds.
loading channel 253 / 384
Elapsed time is 0.034124 seconds.
loading channel 254 / 384
Elapsed time is 0.030693 seconds.
loading channel 255 / 384
Elapsed time is 0.029976 seconds.
loading channel 256 / 384
Elapsed time is 0.024889 seconds.
>> allad

allad =

  1x384 cell array

  Columns 1 through 4

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 5 through 8

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 9 through 12

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 13 through 16

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 17 through 20

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 21 through 24

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 25 through 28

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 29 through 32

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 33 through 36

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 37 through 40

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 41 through 44

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 45 through 48

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 49 through 52

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 53 through 56

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 57 through 60

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 61 through 64

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 65 through 68

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 69 through 72

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 73 through 76

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 77 through 80

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 81 through 84

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 85 through 88

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 89 through 92

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 93 through 96

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 97 through 100

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 101 through 104

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 105 through 108

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 109 through 112

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 113 through 116

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 117 through 120

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 121 through 124

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 125 through 128

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 129 through 131

    {1000000x1 double}    {1000000x1 double}    {1000000x1 double}

  Columns 132 through 134

    {1000000x1 double}    {1000000x1 double}    {1000000x1 double}

  Columns 135 through 137

    {1000000x1 double}    {1000000x1 double}    {1000000x1 double}

  Columns 138 through 140

    {1000000x1 double}    {1000000x1 double}    {1000000x1 double}

  Columns 141 through 143

    {1000000x1 double}    {1000000x1 double}    {1000000x1 double}

  Columns 144 through 146

    {1000000x1 double}    {1000000x1 double}    {1000000x1 double}

  Columns 147 through 149

    {1000000x1 double}    {1000000x1 double}    {1000000x1 double}

  Columns 150 through 152

    {1000000x1 double}    {1000000x1 double}    {1000000x1 double}

  Columns 153 through 155

    {1000000x1 double}    {1000000x1 double}    {1000000x1 double}

  Columns 156 through 158

    {1000000x1 double}    {1000000x1 double}    {1000000x1 double}

  Columns 159 through 161

    {1000000x1 double}    {1000000x1 double}    {1000000x1 double}

  Columns 162 through 164

    {1000000x1 double}    {1000000x1 double}    {1000000x1 double}

  Columns 165 through 167

    {1000000x1 double}    {1000000x1 double}    {1000000x1 double}

  Columns 168 through 170

    {1000000x1 double}    {1000000x1 double}    {1000000x1 double}

  Columns 171 through 173

    {1000000x1 double}    {1000000x1 double}    {1000000x1 double}

  Columns 174 through 176

    {1000000x1 double}    {1000000x1 double}    {1000000x1 double}

  Columns 177 through 179

    {1000000x1 double}    {1000000x1 double}    {1000000x1 double}

  Columns 180 through 182

    {1000000x1 double}    {1000000x1 double}    {1000000x1 double}

  Columns 183 through 185

    {1000000x1 double}    {1000000x1 double}    {1000000x1 double}

  Columns 186 through 188

    {1000000x1 double}    {1000000x1 double}    {1000000x1 double}

  Columns 189 through 191

    {1000000x1 double}    {1000000x1 double}    {1000000x1 double}

  Columns 192 through 194

    {1000000x1 double}    {1000000x1 double}    {1000000x1 double}

  Columns 195 through 197

    {1000000x1 double}    {1000000x1 double}    {1000000x1 double}

  Columns 198 through 200

    {1000000x1 double}    {1000000x1 double}    {1000000x1 double}

  Columns 201 through 203

    {1000000x1 double}    {1000000x1 double}    {1000000x1 double}

  Columns 204 through 206

    {1000000x1 double}    {1000000x1 double}    {1000000x1 double}

  Columns 207 through 209

    {1000000x1 double}    {1000000x1 double}    {1000000x1 double}

  Columns 210 through 212

    {1000000x1 double}    {1000000x1 double}    {1000000x1 double}

  Columns 213 through 215

    {1000000x1 double}    {1000000x1 double}    {1000000x1 double}

  Columns 216 through 218

    {1000000x1 double}    {1000000x1 double}    {1000000x1 double}

  Columns 219 through 221

    {1000000x1 double}    {1000000x1 double}    {1000000x1 double}

  Columns 222 through 224

    {1000000x1 double}    {1000000x1 double}    {1000000x1 double}

  Columns 225 through 227

    {1000000x1 double}    {1000000x1 double}    {1000000x1 double}

  Columns 228 through 230

    {1000000x1 double}    {1000000x1 double}    {1000000x1 double}

  Columns 231 through 233

    {1000000x1 double}    {1000000x1 double}    {1000000x1 double}

  Columns 234 through 236

    {1000000x1 double}    {1000000x1 double}    {1000000x1 double}

  Columns 237 through 239

    {1000000x1 double}    {1000000x1 double}    {1000000x1 double}

  Columns 240 through 242

    {1000000x1 double}    {1000000x1 double}    {1000000x1 double}

  Columns 243 through 245

    {1000000x1 double}    {1000000x1 double}    {1000000x1 double}

  Columns 246 through 248

    {1000000x1 double}    {1000000x1 double}    {1000000x1 double}

  Columns 249 through 251

    {1000000x1 double}    {1000000x1 double}    {1000000x1 double}

  Columns 252 through 254

    {1000000x1 double}    {1000000x1 double}    {1000000x1 double}

  Columns 255 through 258

    {1000000x1 double}    {1000000x1 double}    {0x0 double}    {0x0 double}

  Columns 259 through 262

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 263 through 266

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 267 through 270

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 271 through 274

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 275 through 278

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 279 through 282

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 283 through 286

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 287 through 290

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 291 through 294

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 295 through 298

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 299 through 302

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 303 through 306

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 307 through 310

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 311 through 314

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 315 through 318

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 319 through 322

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 323 through 326

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 327 through 330

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 331 through 334

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 335 through 338

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 339 through 342

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 343 through 346

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 347 through 350

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 351 through 354

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 355 through 358

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 359 through 362

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 363 through 366

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 367 through 370

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 371 through 374

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 375 through 378

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 379 through 382

    {0x0 double}    {0x0 double}    {0x0 double}    {0x0 double}

  Columns 383 through 384

    {0x0 double}    {0x0 double}

>> 
>> loaded = find(cellfun(@numel, allad, 'uniformoutput',true));
>> loaded

loaded =

  Columns 1 through 13

   129   130   131   132   133   134   135   136   137   138   139   140   141

  Columns 14 through 26

   142   143   144   145   146   147   148   149   150   151   152   153   154

  Columns 27 through 39

   155   156   157   158   159   160   161   162   163   164   165   166   167

  Columns 40 through 52

   168   169   170   171   172   173   174   175   176   177   178   179   180

  Columns 53 through 65

   181   182   183   184   185   186   187   188   189   190   191   192   193

  Columns 66 through 78

   194   195   196   197   198   199   200   201   202   203   204   205   206

  Columns 79 through 91

   207   208   209   210   211   212   213   214   215   216   217   218   219

  Columns 92 through 104

   220   221   222   223   224   225   226   227   228   229   230   231   232

  Columns 105 through 117

   233   234   235   236   237   238   239   240   241   242   243   244   245

  Columns 118 through 128

   246   247   248   249   250   251   252   253   254   255   256

>> 
   129   130   131   132   133   134   135   136   137   138   139   140   141

  Columns 14 through 26

   142   143   144   145   146   147   148   149   150   151   152   153   154

  Columns 27 through 39

   155   156   157   158   159   160   161   162   163   164   165   166   167

  Columns 40 through 52

   168   169   170   171   172   173   174   175   176   177   178   179   180

  Columns 53 through 65

   181   182   183   184   185   186   187   188   189   190   191   192   193

  Columns 66 through 78

   194   195   196   197   198   199   200   201   202   203   204   205   206

  Columns 79 through 91

   207   208   209   210   211   212   213   214   215   216   217   218   219

  Columns 92 through 104

   220   221   222   223   224   225   226   227   228   229   230   231   232

  Columns 105 through 117

   233   234   235   236   237   238   239   240   241   242   243   244   245

  Columns 118 through 128

   246   247   248   249   250   251   252   253   254   255   256

>> 

>>    129   130   131   132   133   134   135   136   137   138   139   140   141
    129   130   131   132   133   134   135   136   137   138   139   140   141
          |
Error: Unexpected MATLAB expression.
 
>> 
>>   Columns 14 through 26
Undefined function or variable 'Columns'.
 
>> 
>>    142   143   144   145   146   147   148   149   150   151   152   153   154
    142   143   144   145   146   147   148   149   150   151   152   153   154
          |
Error: Unexpected MATLAB expression.
 
>> 
>>   Columns 27 through 39
Undefined function or variable 'Columns'.
 
  C-c C-c>> 
>> 
>> 
>> 
>> loaded = find(cellfun(@numel, allad, 'uniformoutput',true));
>> loaded = find(cellfun(@numel, allad));
>> 
>> 
>> tmp = horzcat(allad {loaded} );
>> size(tmp)

ans =

     1000000         128

>> plot(tmp(1:1e5,:))
>> clf
>> plot(tmp(1:1e5,:))
>> 