function Read_thermistor(path,filename,status)

%-------------------------------------------------------------------------
% Instructions:
%
% At the Matlab command prompt, 
% set variable called path to current directory, e..g.
%  >> path = '/data/07study/temperature/LHP/LHP002/'
%
% set a variable called filename to the name of the text file which has the data, e.g.
%  >> filename = 'LHP002A000T001.csv'
%
% status = 0 or 1, '0' for temp only, '1' for temp and pressure
%
% Then call this function:
%  >> temptime(filename)
%
%The program will read data, convert all fields to numbers, and save file
%as:  
%    'OAHUAXA015T001.MAT'
%
%==========================================================================

%READ DATA INTO FILES------------------------------------------------------

if status==0
[temperature,date,time] = textread([path,'/',filename],'%f%s%s','delimiter',',','headerlines',37);
elseif status==1
[temperature,pressure,date,time] = textread([path,'/',filename],'%f%f%s%s','delimiter',',','headerlines',46);
A.pressure=pressure;
end

%Fill in missing data with NaN

ff = find(temperature == 9999);
temperature(ff) = NaN;

%Convert date string to number in julian day and mtime

date2 = char(date);
time2 = char(time);

hh = str2num(time2(:,1:2));
mm = str2num(time2(:,4:5));
ss = str2num(time2(:,7:8));

mtime = datenum(date)+datenum(0,0,0,hh,mm,ss);
julday=mtime-datenum(2009,0,0);

sitecode = filename(1:6);
file_name = [filename(1:end-2),'TXT'];

%Set site-specific fields

A.sitecode = sitecode;
A.description = input('enter site title: ','s');
%A.lat = input('enter lat: ','s');
%A.lon = input('enter long: ','s');
A.waterdepth = input('enter depth: '); 
A.depth_logger = input('enter thermistor depth: ');
A.height = A.waterdepth - A.depth_logger;
A.filename = [filename(1:end-3)];
A.temperature = temperature;

A.mtime = mtime;
A.julday=julday;
A.name = 'temp_qc';

A = orderfields(A);
temp_data = A;


eval(['save ',path,'/',A.filename,'MAT',' temp_data']);
    
end







