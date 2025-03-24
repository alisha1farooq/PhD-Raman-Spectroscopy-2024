% WiREDataUnit  Enumeration used by WiRE to indicate units of measurement

% Copyright Notice:
% (c) 2012 - 2013 Renishaw plc. All rights reserved.

classdef WiREDataUnit < uint32
    enumeration
        Arbitrary          (0),
        RamanShift         (1),
        Wavenumber         (2),
        Nanometre          (3),
        ElectronVolt       (4),
        Micron             (5),
        Counts             (6),
        Electrons          (7),
        Millimetres        (8),
        Metres             (9),
        Kelvin            (10),
        Pascal            (11),
        Seconds           (12),
        Milliseconds      (13),
        Hours             (14),
        Days              (15),
        Pixels            (16),
        Intensity         (17),
        RelativeIntensity (18),
        Degrees           (19),
        Radians           (20),
        Celsius           (21),
        Fahrenheit        (22),
        KelvinPerMinute   (23),
        FileTime          (24)
    end;
end