pathname='C:\Dokumente und Einstellungen\Besitzer\Desktop\gbpusd data\';
filename='GBPUSD_2001_2006.txt';
%<TICKER>,<DTYYYYMMDD>,<TIME>,<OPEN>,<HIGH>,<LOW>,<CLOSE>
[ticker,bigdate,itime,po,ph,pl,pc] = textread([pathname filename],'%6c,%8n,%6n,%f,%f,%f,%f');
