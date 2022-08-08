EESchema Schematic File Version 4
EELAYER 30 0
EELAYER END
$Descr A4 11693 8268
encoding utf-8
Sheet 1 1
Title ""
Date ""
Rev ""
Comp ""
Comment1 ""
Comment2 ""
Comment3 ""
Comment4 ""
$EndDescr
$Comp
L 74xx:74HC4051 U1
U 1 1 622D14F8
P 2600 2900
F 0 "U1" H 2807 3581 50  0000 C CNN
F 1 "74HC4051" H 2886 3490 50  0000 C CNN
F 2 "Package_SO:SOP-16_4.55x10.3mm_P1.27mm" H 2600 2500 50  0001 C CNN
F 3 "http://www.ti.com/lit/ds/symlink/cd74hc4051.pdf" H 2600 2500 50  0001 C CNN
	1    2600 2900
	1    0    0    -1  
$EndComp
$Comp
L 74xx:74HC4051 U2
U 1 1 622FAE29
P 8200 2900
F 0 "U2" H 8486 3581 50  0000 C CNN
F 1 "74HC4051" H 8487 3490 50  0000 C CNN
F 2 "Package_SO:SOP-16_4.55x10.3mm_P1.27mm" H 8200 2500 50  0001 C CNN
F 3 "http://www.ti.com/lit/ds/symlink/cd74hc4051.pdf" H 8200 2500 50  0001 C CNN
	1    8200 2900
	-1   0    0    -1  
$EndComp
$Comp
L Device:C C28
U 1 1 622FAE2F
P 9900 2950
F 0 "C28" H 10015 2996 50  0000 L CNN
F 1 "100n" H 10015 2905 50  0000 L CNN
F 2 "Capacitor_THT:C_Rect_L7.0mm_W2.0mm_P5.00mm" H 9938 2800 50  0001 C CNN
F 3 "~" H 9900 2950 50  0001 C CNN
	1    9900 2950
	-1   0    0    1   
$EndComp
$Comp
L Device:C C26
U 1 1 622FAE35
P 9800 1850
F 0 "C26" H 9915 1896 50  0000 L CNN
F 1 "100n" H 9915 1805 50  0000 L CNN
F 2 "Capacitor_THT:C_Rect_L7.0mm_W2.0mm_P5.00mm" H 9838 1700 50  0001 C CNN
F 3 "~" H 9800 1850 50  0001 C CNN
	1    9800 1850
	0    -1   -1   0   
$EndComp
$Comp
L Device:C C27
U 1 1 622FAE3B
P 9900 2400
F 0 "C27" H 10015 2446 50  0000 L CNN
F 1 "100n" H 10015 2355 50  0000 L CNN
F 2 "Capacitor_THT:C_Rect_L7.0mm_W2.0mm_P5.00mm" H 9938 2250 50  0001 C CNN
F 3 "~" H 9900 2400 50  0001 C CNN
	1    9900 2400
	-1   0    0    1   
$EndComp
$Comp
L Device:R R7
U 1 1 622FAE41
P 9550 2900
F 0 "R7" H 9620 2946 50  0000 L CNN
F 1 "15k" H 9620 2855 50  0000 L CNN
F 2 "Resistor_THT:R_Axial_DIN0207_L6.3mm_D2.5mm_P2.54mm_Vertical" V 9480 2900 50  0001 C CNN
F 3 "~" H 9550 2900 50  0001 C CNN
	1    9550 2900
	-1   0    0    1   
$EndComp
$Comp
L Device:R R8
U 1 1 622FAE47
P 9550 2400
F 0 "R8" H 9620 2446 50  0000 L CNN
F 1 "15k" H 9620 2355 50  0000 L CNN
F 2 "Resistor_THT:R_Axial_DIN0207_L6.3mm_D2.5mm_P2.54mm_Vertical" V 9480 2400 50  0001 C CNN
F 3 "~" H 9550 2400 50  0001 C CNN
	1    9550 2400
	-1   0    0    1   
$EndComp
$Comp
L Device:R R6
U 1 1 622FAE4D
P 9150 2600
F 0 "R6" H 9220 2646 50  0000 L CNN
F 1 "10k" H 9220 2555 50  0000 L CNN
F 2 "Resistor_THT:R_Axial_DIN0207_L6.3mm_D2.5mm_P2.54mm_Vertical" V 9080 2600 50  0001 C CNN
F 3 "~" H 9150 2600 50  0001 C CNN
	1    9150 2600
	0    1    1    0   
$EndComp
$Comp
L Device:C C4
U 1 1 62327C88
P 4250 1400
F 0 "C4" H 4365 1446 50  0000 L CNN
F 1 "100n" H 4365 1355 50  0000 L CNN
F 2 "Capacitor_THT:C_Rect_L7.0mm_W2.0mm_P5.00mm" H 4288 1250 50  0001 C CNN
F 3 "~" H 4250 1400 50  0001 C CNN
	1    4250 1400
	0    -1   -1   0   
$EndComp
$Comp
L Device:C C6
U 1 1 62329CB9
P 4800 1550
F 0 "C6" H 4915 1596 50  0000 L CNN
F 1 "1500" H 4915 1505 50  0000 L CNN
F 2 "Capacitor_THT:C_Disc_D4.3mm_W1.9mm_P5.00mm" H 4838 1400 50  0001 C CNN
F 3 "~" H 4800 1550 50  0001 C CNN
	1    4800 1550
	1    0    0    -1  
$EndComp
$Comp
L Device:C C7
U 1 1 6232A9A1
P 6100 1550
F 0 "C7" H 6215 1596 50  0000 L CNN
F 1 "1500" H 6215 1505 50  0000 L CNN
F 2 "Capacitor_THT:C_Disc_D4.3mm_W1.9mm_P5.00mm" H 6138 1400 50  0001 C CNN
F 3 "~" H 6100 1550 50  0001 C CNN
	1    6100 1550
	1    0    0    -1  
$EndComp
$Comp
L Device:C C8
U 1 1 6232AEBB
P 6450 1400
F 0 "C8" H 6565 1446 50  0000 L CNN
F 1 "100n" H 6565 1355 50  0000 L CNN
F 2 "Capacitor_THT:C_Rect_L7.0mm_W2.0mm_P5.00mm" H 6488 1250 50  0001 C CNN
F 3 "~" H 6450 1400 50  0001 C CNN
	1    6450 1400
	0    -1   -1   0   
$EndComp
$Comp
L Device:C C5
U 1 1 6232EC8E
P 5150 1400
F 0 "C5" H 5265 1446 50  0000 L CNN
F 1 "680" H 5265 1355 50  0000 L CNN
F 2 "Capacitor_THT:C_Disc_D4.3mm_W1.9mm_P5.00mm" H 5188 1250 50  0001 C CNN
F 3 "~" H 5150 1400 50  0001 C CNN
	1    5150 1400
	0    -1   -1   0   
$EndComp
$Comp
L Device:L L1
U 1 1 623348B0
P 4550 1550
F 0 "L1" H 4485 1596 50  0000 R CNN
F 1 "2.2uH" H 4485 1505 50  0000 R CNN
F 2 "Inductor_THT:L_Axial_L5.3mm_D2.2mm_P2.54mm_Vertical_Vishay_IM-1" H 4550 1550 50  0001 C CNN
F 3 "~" H 4550 1550 50  0001 C CNN
	1    4550 1550
	1    0    0    -1  
$EndComp
$Comp
L Device:L L3
U 1 1 62334E2A
P 5850 1550
F 0 "L3" H 5785 1596 50  0000 R CNN
F 1 "2.2uH" H 5785 1505 50  0000 R CNN
F 2 "Inductor_THT:L_Axial_L5.3mm_D2.2mm_P2.54mm_Vertical_Vishay_IM-1" H 5850 1550 50  0001 R CNN
F 3 "~" H 5850 1550 50  0001 C CNN
	1    5850 1550
	1    0    0    -1  
$EndComp
$Comp
L Device:L L2
U 1 1 62335445
P 5600 1400
F 0 "L2" H 5653 1446 50  0000 L CNN
F 1 "4.7uH" H 5653 1355 50  0000 L CNN
F 2 "Inductor_THT:L_Axial_L5.3mm_D2.2mm_P2.54mm_Vertical_Vishay_IM-1" H 5600 1400 50  0001 C CNN
F 3 "~" H 5600 1400 50  0001 C CNN
	1    5600 1400
	0    -1   -1   0   
$EndComp
Wire Wire Line
	5750 1400 5850 1400
Wire Wire Line
	6100 1400 6300 1400
Wire Wire Line
	6100 1400 5850 1400
Connection ~ 6100 1400
Connection ~ 5850 1400
Wire Wire Line
	5450 1400 5300 1400
Wire Wire Line
	5000 1400 4800 1400
Wire Wire Line
	4800 1400 4550 1400
Connection ~ 4800 1400
Wire Wire Line
	4550 1400 4400 1400
Connection ~ 4550 1400
$Comp
L power:GNDREF #PWR09
U 1 1 622EEAEB
P 4550 1700
F 0 "#PWR09" H 4550 1450 50  0001 C CNN
F 1 "GNDREF" H 4555 1527 50  0001 C CNN
F 2 "" H 4550 1700 50  0001 C CNN
F 3 "" H 4550 1700 50  0001 C CNN
	1    4550 1700
	1    0    0    -1  
$EndComp
$Comp
L power:GNDREF #PWR013
U 1 1 622EF224
P 4800 1700
F 0 "#PWR013" H 4800 1450 50  0001 C CNN
F 1 "GNDREF" H 4805 1527 50  0001 C CNN
F 2 "" H 4800 1700 50  0001 C CNN
F 3 "" H 4800 1700 50  0001 C CNN
	1    4800 1700
	1    0    0    -1  
$EndComp
$Comp
L power:GNDREF #PWR017
U 1 1 622EF644
P 5850 1700
F 0 "#PWR017" H 5850 1450 50  0001 C CNN
F 1 "GNDREF" H 5855 1527 50  0001 C CNN
F 2 "" H 5850 1700 50  0001 C CNN
F 3 "" H 5850 1700 50  0001 C CNN
	1    5850 1700
	1    0    0    -1  
$EndComp
$Comp
L power:GNDREF #PWR021
U 1 1 622F1304
P 6100 1700
F 0 "#PWR021" H 6100 1450 50  0001 C CNN
F 1 "GNDREF" H 6105 1527 50  0001 C CNN
F 2 "" H 6100 1700 50  0001 C CNN
F 3 "" H 6100 1700 50  0001 C CNN
	1    6100 1700
	1    0    0    -1  
$EndComp
Text Label 5000 1100 0    50   ~ 0
160-80m
$Comp
L Device:C C9
U 1 1 622F8E59
P 4250 2400
F 0 "C9" H 4365 2446 50  0000 L CNN
F 1 "100n" H 4365 2355 50  0000 L CNN
F 2 "Capacitor_THT:C_Rect_L7.0mm_W2.0mm_P5.00mm" H 4288 2250 50  0001 C CNN
F 3 "~" H 4250 2400 50  0001 C CNN
	1    4250 2400
	0    -1   -1   0   
$EndComp
$Comp
L Device:C C10
U 1 1 622F8E5F
P 4800 2550
F 0 "C10" H 4915 2596 50  0000 L CNN
F 1 "820" H 4915 2505 50  0000 L CNN
F 2 "Capacitor_THT:C_Disc_D4.3mm_W1.9mm_P5.00mm" H 4838 2400 50  0001 C CNN
F 3 "~" H 4800 2550 50  0001 C CNN
	1    4800 2550
	1    0    0    -1  
$EndComp
$Comp
L Device:C C12
U 1 1 622F8E65
P 6100 2550
F 0 "C12" H 6215 2596 50  0000 L CNN
F 1 "820" H 6215 2505 50  0000 L CNN
F 2 "Capacitor_THT:C_Disc_D4.3mm_W1.9mm_P5.00mm" H 6138 2400 50  0001 C CNN
F 3 "~" H 6100 2550 50  0001 C CNN
	1    6100 2550
	1    0    0    -1  
$EndComp
$Comp
L Device:C C13
U 1 1 622F8E6B
P 6450 2400
F 0 "C13" H 6565 2446 50  0000 L CNN
F 1 "100n" H 6565 2355 50  0000 L CNN
F 2 "Capacitor_THT:C_Rect_L7.0mm_W2.0mm_P5.00mm" H 6488 2250 50  0001 C CNN
F 3 "~" H 6450 2400 50  0001 C CNN
	1    6450 2400
	0    -1   -1   0   
$EndComp
$Comp
L Device:C C11
U 1 1 622F8E71
P 5150 2400
F 0 "C11" H 5265 2446 50  0000 L CNN
F 1 "390" H 5265 2355 50  0000 L CNN
F 2 "Capacitor_THT:C_Disc_D4.3mm_W1.9mm_P5.00mm" H 5188 2250 50  0001 C CNN
F 3 "~" H 5150 2400 50  0001 C CNN
	1    5150 2400
	0    -1   -1   0   
$EndComp
$Comp
L Device:L L4
U 1 1 622F8E77
P 4550 2550
F 0 "L4" H 4485 2596 50  0000 R CNN
F 1 "1uH" H 4485 2505 50  0000 R CNN
F 2 "Inductor_THT:L_Axial_L5.3mm_D2.2mm_P2.54mm_Vertical_Vishay_IM-1" H 4550 2550 50  0001 C CNN
F 3 "~" H 4550 2550 50  0001 C CNN
	1    4550 2550
	1    0    0    -1  
$EndComp
$Comp
L Device:L L6
U 1 1 622F8E7D
P 5850 2550
F 0 "L6" H 5785 2596 50  0000 R CNN
F 1 "1uH" H 5785 2505 50  0000 R CNN
F 2 "Inductor_THT:L_Axial_L5.3mm_D2.2mm_P2.54mm_Vertical_Vishay_IM-1" H 5850 2550 50  0001 R CNN
F 3 "~" H 5850 2550 50  0001 C CNN
	1    5850 2550
	1    0    0    -1  
$EndComp
$Comp
L Device:L L5
U 1 1 622F8E83
P 5600 2400
F 0 "L5" H 5653 2446 50  0000 L CNN
F 1 "2.2uH" H 5653 2355 50  0000 L CNN
F 2 "Inductor_THT:L_Axial_L5.3mm_D2.2mm_P2.54mm_Vertical_Vishay_IM-1" H 5600 2400 50  0001 C CNN
F 3 "~" H 5600 2400 50  0001 C CNN
	1    5600 2400
	0    -1   -1   0   
$EndComp
Wire Wire Line
	5750 2400 5850 2400
Wire Wire Line
	6100 2400 6300 2400
Wire Wire Line
	6100 2400 5850 2400
Connection ~ 6100 2400
Connection ~ 5850 2400
Wire Wire Line
	5450 2400 5300 2400
Wire Wire Line
	5000 2400 4800 2400
Wire Wire Line
	4800 2400 4550 2400
Connection ~ 4800 2400
Wire Wire Line
	4550 2400 4400 2400
Connection ~ 4550 2400
$Comp
L power:GNDREF #PWR010
U 1 1 622F8E94
P 4550 2700
F 0 "#PWR010" H 4550 2450 50  0001 C CNN
F 1 "GNDREF" H 4555 2527 50  0001 C CNN
F 2 "" H 4550 2700 50  0001 C CNN
F 3 "" H 4550 2700 50  0001 C CNN
	1    4550 2700
	1    0    0    -1  
$EndComp
$Comp
L power:GNDREF #PWR014
U 1 1 622F8E9A
P 4800 2700
F 0 "#PWR014" H 4800 2450 50  0001 C CNN
F 1 "GNDREF" H 4805 2527 50  0001 C CNN
F 2 "" H 4800 2700 50  0001 C CNN
F 3 "" H 4800 2700 50  0001 C CNN
	1    4800 2700
	1    0    0    -1  
$EndComp
$Comp
L power:GNDREF #PWR018
U 1 1 622F8EA0
P 5850 2700
F 0 "#PWR018" H 5850 2450 50  0001 C CNN
F 1 "GNDREF" H 5855 2527 50  0001 C CNN
F 2 "" H 5850 2700 50  0001 C CNN
F 3 "" H 5850 2700 50  0001 C CNN
	1    5850 2700
	1    0    0    -1  
$EndComp
$Comp
L power:GNDREF #PWR022
U 1 1 622F8EA6
P 6100 2700
F 0 "#PWR022" H 6100 2450 50  0001 C CNN
F 1 "GNDREF" H 6105 2527 50  0001 C CNN
F 2 "" H 6100 2700 50  0001 C CNN
F 3 "" H 6100 2700 50  0001 C CNN
	1    6100 2700
	1    0    0    -1  
$EndComp
Text Label 5000 2100 0    50   ~ 0
40m
$Comp
L Device:C C14
U 1 1 622FA8B1
P 4200 3450
F 0 "C14" H 4315 3496 50  0000 L CNN
F 1 "100n" H 4315 3405 50  0000 L CNN
F 2 "Capacitor_THT:C_Rect_L7.0mm_W2.0mm_P5.00mm" H 4238 3300 50  0001 C CNN
F 3 "~" H 4200 3450 50  0001 C CNN
	1    4200 3450
	0    -1   -1   0   
$EndComp
$Comp
L Device:C C15
U 1 1 622FA8B7
P 4750 3600
F 0 "C15" H 4865 3646 50  0000 L CNN
F 1 "390" H 4865 3555 50  0000 L CNN
F 2 "Capacitor_THT:C_Disc_D4.3mm_W1.9mm_P5.00mm" H 4788 3450 50  0001 C CNN
F 3 "~" H 4750 3600 50  0001 C CNN
	1    4750 3600
	1    0    0    -1  
$EndComp
$Comp
L Device:C C17
U 1 1 622FA8BD
P 6050 3600
F 0 "C17" H 6165 3646 50  0000 L CNN
F 1 "390" H 6165 3555 50  0000 L CNN
F 2 "Capacitor_THT:C_Disc_D4.3mm_W1.9mm_P5.00mm" H 6088 3450 50  0001 C CNN
F 3 "~" H 6050 3600 50  0001 C CNN
	1    6050 3600
	1    0    0    -1  
$EndComp
$Comp
L Device:C C18
U 1 1 622FA8C3
P 6400 3450
F 0 "C18" H 6515 3496 50  0000 L CNN
F 1 "100n" H 6515 3405 50  0000 L CNN
F 2 "Capacitor_THT:C_Rect_L7.0mm_W2.0mm_P5.00mm" H 6438 3300 50  0001 C CNN
F 3 "~" H 6400 3450 50  0001 C CNN
	1    6400 3450
	0    -1   -1   0   
$EndComp
$Comp
L Device:C C16
U 1 1 622FA8C9
P 5100 3450
F 0 "C16" H 5215 3496 50  0000 L CNN
F 1 "180" H 5215 3405 50  0000 L CNN
F 2 "Capacitor_THT:C_Disc_D4.3mm_W1.9mm_P5.00mm" H 5138 3300 50  0001 C CNN
F 3 "~" H 5100 3450 50  0001 C CNN
	1    5100 3450
	0    -1   -1   0   
$EndComp
$Comp
L Device:L L7
U 1 1 622FA8CF
P 4500 3600
F 0 "L7" H 4435 3646 50  0000 R CNN
F 1 "0.47uH" H 4435 3555 50  0000 R CNN
F 2 "Inductor_THT:L_Axial_L5.3mm_D2.2mm_P2.54mm_Vertical_Vishay_IM-1" H 4500 3600 50  0001 C CNN
F 3 "~" H 4500 3600 50  0001 C CNN
	1    4500 3600
	1    0    0    -1  
$EndComp
$Comp
L Device:L L9
U 1 1 622FA8D5
P 5800 3600
F 0 "L9" H 5735 3646 50  0000 R CNN
F 1 "0.47uH" H 5735 3555 50  0000 R CNN
F 2 "Inductor_THT:L_Axial_L5.3mm_D2.2mm_P2.54mm_Vertical_Vishay_IM-1" H 5800 3600 50  0001 R CNN
F 3 "~" H 5800 3600 50  0001 C CNN
	1    5800 3600
	1    0    0    -1  
$EndComp
$Comp
L Device:L L8
U 1 1 622FA8DB
P 5550 3450
F 0 "L8" H 5603 3496 50  0000 L CNN
F 1 "1uH" H 5603 3405 50  0000 L CNN
F 2 "Inductor_THT:L_Axial_L5.3mm_D2.2mm_P2.54mm_Vertical_Vishay_IM-1" H 5550 3450 50  0001 C CNN
F 3 "~" H 5550 3450 50  0001 C CNN
	1    5550 3450
	0    -1   -1   0   
$EndComp
Wire Wire Line
	5700 3450 5800 3450
Wire Wire Line
	6050 3450 6250 3450
Wire Wire Line
	6050 3450 5800 3450
Connection ~ 6050 3450
Connection ~ 5800 3450
Wire Wire Line
	5400 3450 5250 3450
Wire Wire Line
	4950 3450 4750 3450
Wire Wire Line
	4750 3450 4500 3450
Connection ~ 4750 3450
Wire Wire Line
	4500 3450 4350 3450
Connection ~ 4500 3450
$Comp
L power:GNDREF #PWR07
U 1 1 622FA8EC
P 4500 3750
F 0 "#PWR07" H 4500 3500 50  0001 C CNN
F 1 "GNDREF" H 4505 3577 50  0001 C CNN
F 2 "" H 4500 3750 50  0001 C CNN
F 3 "" H 4500 3750 50  0001 C CNN
	1    4500 3750
	1    0    0    -1  
$EndComp
$Comp
L power:GNDREF #PWR011
U 1 1 622FA8F2
P 4750 3750
F 0 "#PWR011" H 4750 3500 50  0001 C CNN
F 1 "GNDREF" H 4755 3577 50  0001 C CNN
F 2 "" H 4750 3750 50  0001 C CNN
F 3 "" H 4750 3750 50  0001 C CNN
	1    4750 3750
	1    0    0    -1  
$EndComp
$Comp
L power:GNDREF #PWR015
U 1 1 622FA8F8
P 5800 3750
F 0 "#PWR015" H 5800 3500 50  0001 C CNN
F 1 "GNDREF" H 5805 3577 50  0001 C CNN
F 2 "" H 5800 3750 50  0001 C CNN
F 3 "" H 5800 3750 50  0001 C CNN
	1    5800 3750
	1    0    0    -1  
$EndComp
$Comp
L power:GNDREF #PWR019
U 1 1 622FA8FE
P 6050 3750
F 0 "#PWR019" H 6050 3500 50  0001 C CNN
F 1 "GNDREF" H 6055 3577 50  0001 C CNN
F 2 "" H 6050 3750 50  0001 C CNN
F 3 "" H 6050 3750 50  0001 C CNN
	1    6050 3750
	1    0    0    -1  
$EndComp
Text Label 4950 3150 0    50   ~ 0
30-20m
$Comp
L Device:C C19
U 1 1 622FFBE4
P 4200 4500
F 0 "C19" H 4315 4546 50  0000 L CNN
F 1 "100n" H 4315 4455 50  0000 L CNN
F 2 "Capacitor_THT:C_Rect_L7.0mm_W2.0mm_P5.00mm" H 4238 4350 50  0001 C CNN
F 3 "~" H 4200 4500 50  0001 C CNN
	1    4200 4500
	0    -1   -1   0   
$EndComp
$Comp
L Device:C C20
U 1 1 622FFBEA
P 4750 4650
F 0 "C20" H 4865 4696 50  0000 L CNN
F 1 "180" H 4865 4605 50  0000 L CNN
F 2 "Capacitor_THT:C_Disc_D4.3mm_W1.9mm_P5.00mm" H 4788 4500 50  0001 C CNN
F 3 "~" H 4750 4650 50  0001 C CNN
	1    4750 4650
	1    0    0    -1  
$EndComp
$Comp
L Device:C C22
U 1 1 622FFBF0
P 6050 4650
F 0 "C22" H 6165 4696 50  0000 L CNN
F 1 "180" H 6165 4605 50  0000 L CNN
F 2 "Capacitor_THT:C_Disc_D4.3mm_W1.9mm_P5.00mm" H 6088 4500 50  0001 C CNN
F 3 "~" H 6050 4650 50  0001 C CNN
	1    6050 4650
	1    0    0    -1  
$EndComp
$Comp
L Device:C C23
U 1 1 622FFBF6
P 6400 4500
F 0 "C23" H 6515 4546 50  0000 L CNN
F 1 "100n" H 6515 4455 50  0000 L CNN
F 2 "Capacitor_THT:C_Rect_L7.0mm_W2.0mm_P5.00mm" H 6438 4350 50  0001 C CNN
F 3 "~" H 6400 4500 50  0001 C CNN
	1    6400 4500
	0    -1   -1   0   
$EndComp
$Comp
L Device:C C21
U 1 1 622FFBFC
P 5100 4500
F 0 "C21" H 5215 4546 50  0000 L CNN
F 1 "82" H 5215 4455 50  0000 L CNN
F 2 "Capacitor_THT:C_Disc_D4.3mm_W1.9mm_P5.00mm" H 5138 4350 50  0001 C CNN
F 3 "~" H 5100 4500 50  0001 C CNN
	1    5100 4500
	0    -1   -1   0   
$EndComp
$Comp
L Device:L L10
U 1 1 622FFC02
P 4500 4650
F 0 "L10" H 4435 4696 50  0000 R CNN
F 1 "0.33uH" H 4435 4605 50  0000 R CNN
F 2 "Inductor_THT:L_Axial_L5.3mm_D2.2mm_P2.54mm_Vertical_Vishay_IM-1" H 4500 4650 50  0001 C CNN
F 3 "~" H 4500 4650 50  0001 C CNN
	1    4500 4650
	1    0    0    -1  
$EndComp
$Comp
L Device:L L12
U 1 1 622FFC08
P 5800 4650
F 0 "L12" H 5735 4696 50  0000 R CNN
F 1 "0.33uH" H 5735 4605 50  0000 R CNN
F 2 "Inductor_THT:L_Axial_L5.3mm_D2.2mm_P2.54mm_Vertical_Vishay_IM-1" H 5800 4650 50  0001 R CNN
F 3 "~" H 5800 4650 50  0001 C CNN
	1    5800 4650
	1    0    0    -1  
$EndComp
$Comp
L Device:L L11
U 1 1 622FFC0E
P 5550 4500
F 0 "L11" H 5603 4546 50  0000 L CNN
F 1 "0.68uH" H 5603 4455 50  0000 L CNN
F 2 "Inductor_THT:L_Axial_L5.3mm_D2.2mm_P2.54mm_Vertical_Vishay_IM-1" H 5550 4500 50  0001 C CNN
F 3 "~" H 5550 4500 50  0001 C CNN
	1    5550 4500
	0    -1   -1   0   
$EndComp
Wire Wire Line
	5700 4500 5800 4500
Wire Wire Line
	6050 4500 6250 4500
Wire Wire Line
	6050 4500 5800 4500
Connection ~ 6050 4500
Connection ~ 5800 4500
Wire Wire Line
	5400 4500 5250 4500
Wire Wire Line
	4950 4500 4750 4500
Wire Wire Line
	4750 4500 4500 4500
Connection ~ 4750 4500
Wire Wire Line
	4500 4500 4350 4500
Connection ~ 4500 4500
$Comp
L power:GNDREF #PWR08
U 1 1 622FFC1F
P 4500 4800
F 0 "#PWR08" H 4500 4550 50  0001 C CNN
F 1 "GNDREF" H 4505 4627 50  0001 C CNN
F 2 "" H 4500 4800 50  0001 C CNN
F 3 "" H 4500 4800 50  0001 C CNN
	1    4500 4800
	1    0    0    -1  
$EndComp
$Comp
L power:GNDREF #PWR012
U 1 1 622FFC25
P 4750 4800
F 0 "#PWR012" H 4750 4550 50  0001 C CNN
F 1 "GNDREF" H 4755 4627 50  0001 C CNN
F 2 "" H 4750 4800 50  0001 C CNN
F 3 "" H 4750 4800 50  0001 C CNN
	1    4750 4800
	1    0    0    -1  
$EndComp
$Comp
L power:GNDREF #PWR016
U 1 1 622FFC2B
P 5800 4800
F 0 "#PWR016" H 5800 4550 50  0001 C CNN
F 1 "GNDREF" H 5805 4627 50  0001 C CNN
F 2 "" H 5800 4800 50  0001 C CNN
F 3 "" H 5800 4800 50  0001 C CNN
	1    5800 4800
	1    0    0    -1  
$EndComp
$Comp
L power:GNDREF #PWR020
U 1 1 622FFC31
P 6050 4800
F 0 "#PWR020" H 6050 4550 50  0001 C CNN
F 1 "GNDREF" H 6055 4627 50  0001 C CNN
F 2 "" H 6050 4800 50  0001 C CNN
F 3 "" H 6050 4800 50  0001 C CNN
	1    6050 4800
	1    0    0    -1  
$EndComp
Text Label 4950 4200 0    50   ~ 0
16-10m
Wire Wire Line
	3000 2600 3250 2600
Wire Wire Line
	3250 2600 3250 1400
Wire Wire Line
	3250 1400 4100 1400
Wire Wire Line
	3000 2700 3400 2700
Wire Wire Line
	3400 2700 3400 2400
Wire Wire Line
	3400 2400 4100 2400
Wire Wire Line
	3000 2800 3400 2800
Wire Wire Line
	3400 2800 3400 3450
Wire Wire Line
	3400 3450 4050 3450
Wire Wire Line
	3000 2900 3300 2900
Wire Wire Line
	3300 2900 3300 4500
Wire Wire Line
	3300 4500 4050 4500
$Comp
L Device:C C24
U 1 1 6236E849
P 4950 5450
F 0 "C24" H 5065 5496 50  0000 L CNN
F 1 "100n" H 5065 5405 50  0000 L CNN
F 2 "Capacitor_THT:C_Rect_L7.0mm_W2.0mm_P5.00mm" H 4988 5300 50  0001 C CNN
F 3 "~" H 4950 5450 50  0001 C CNN
	1    4950 5450
	0    -1   -1   0   
$EndComp
Wire Wire Line
	3000 3000 3200 3000
Wire Wire Line
	3200 5450 4800 5450
Wire Wire Line
	3200 3000 3200 5450
Text Label 4800 5100 0    50   ~ 0
PassThrough
Wire Wire Line
	7800 2600 7450 2600
Wire Wire Line
	7450 2600 7450 1400
Wire Wire Line
	7450 1400 6600 1400
Wire Wire Line
	7800 2700 7350 2700
Wire Wire Line
	7350 2700 7350 2400
Wire Wire Line
	7350 2400 6600 2400
Wire Wire Line
	7800 2800 7350 2800
Wire Wire Line
	7350 2800 7350 3450
Wire Wire Line
	7350 3450 6550 3450
Wire Wire Line
	7800 2900 7450 2900
Wire Wire Line
	7450 2900 7450 4500
Wire Wire Line
	7450 4500 6550 4500
Wire Wire Line
	7800 3000 7550 3000
Wire Wire Line
	7550 3000 7550 5450
Wire Wire Line
	7550 5450 5100 5450
$Comp
L power:GNDREF #PWR023
U 1 1 6239EE98
P 8100 3500
F 0 "#PWR023" H 8100 3250 50  0001 C CNN
F 1 "GNDREF" H 8105 3327 50  0001 C CNN
F 2 "" H 8100 3500 50  0001 C CNN
F 3 "" H 8100 3500 50  0001 C CNN
	1    8100 3500
	1    0    0    -1  
$EndComp
$Comp
L power:GNDREF #PWR024
U 1 1 6239F625
P 8200 3500
F 0 "#PWR024" H 8200 3250 50  0001 C CNN
F 1 "GNDREF" H 8205 3327 50  0001 C CNN
F 2 "" H 8200 3500 50  0001 C CNN
F 3 "" H 8200 3500 50  0001 C CNN
	1    8200 3500
	1    0    0    -1  
$EndComp
$Comp
L Device:R R11
U 1 1 623C1241
P 9150 5800
F 0 "R11" H 9220 5846 50  0000 L CNN
F 1 "1k" H 9220 5755 50  0000 L CNN
F 2 "Resistor_THT:R_Axial_DIN0207_L6.3mm_D2.5mm_P2.54mm_Vertical" V 9080 5800 50  0001 C CNN
F 3 "~" H 9150 5800 50  0001 C CNN
	1    9150 5800
	1    0    0    -1  
$EndComp
$Comp
L Device:R R10
U 1 1 623C19BA
P 8850 5800
F 0 "R10" H 8920 5846 50  0000 L CNN
F 1 "1k" H 8920 5755 50  0000 L CNN
F 2 "Resistor_THT:R_Axial_DIN0207_L6.3mm_D2.5mm_P2.54mm_Vertical" V 8780 5800 50  0001 C CNN
F 3 "~" H 8850 5800 50  0001 C CNN
	1    8850 5800
	1    0    0    -1  
$EndComp
$Comp
L Device:R R9
U 1 1 623C1F1F
P 8550 5800
F 0 "R9" H 8620 5846 50  0000 L CNN
F 1 "1k" H 8620 5755 50  0000 L CNN
F 2 "Resistor_THT:R_Axial_DIN0207_L6.3mm_D2.5mm_P2.54mm_Vertical" V 8480 5800 50  0001 C CNN
F 3 "~" H 8550 5800 50  0001 C CNN
	1    8550 5800
	1    0    0    -1  
$EndComp
Wire Wire Line
	9450 5350 8850 5350
Wire Wire Line
	9150 5650 9150 5150
Connection ~ 9150 5150
Wire Wire Line
	8850 5650 8850 5350
Connection ~ 8850 5350
Wire Wire Line
	8550 5650 8550 5550
Connection ~ 8550 5550
Wire Wire Line
	8550 5950 8850 5950
Connection ~ 8850 5950
$Comp
L power:GNDREF #PWR026
U 1 1 623DA588
P 8850 5950
F 0 "#PWR026" H 8850 5700 50  0001 C CNN
F 1 "GNDREF" H 8855 5777 50  0001 C CNN
F 2 "" H 8850 5950 50  0001 C CNN
F 3 "" H 8850 5950 50  0001 C CNN
	1    8850 5950
	1    0    0    -1  
$EndComp
Wire Wire Line
	8850 5950 9150 5950
Text GLabel 8650 3050 2    50   Input ~ 0
S2
Text GLabel 8650 2900 2    50   Input ~ 0
S1
Text GLabel 8650 2750 2    50   Input ~ 0
S0
Wire Wire Line
	8500 2800 8650 2800
Wire Wire Line
	8650 2800 8650 2750
Wire Wire Line
	8450 2900 8500 2900
Connection ~ 8500 2900
Wire Wire Line
	8500 2900 8650 2900
Wire Wire Line
	8500 3000 8650 3000
Wire Wire Line
	8650 3000 8650 3050
Text GLabel 8250 5150 0    50   Output ~ 0
S0
Text GLabel 8250 5350 0    50   Output ~ 0
S1
Text GLabel 8250 5550 0    50   Output ~ 0
S2
Wire Wire Line
	8250 5150 9150 5150
$Comp
L Switch:SW_Coded_SH-7010 SW1
U 1 1 624264A2
P 10100 5250
F 0 "SW1" H 9770 5296 50  0000 R CNN
F 1 "BCD Switch" H 9770 5205 50  0000 R CNN
F 2 "" H 9800 4800 50  0001 L CNN
F 3 "https://www.nidec-copal-electronics.com/e/catalog/switch/sh-7000.pdf" H 10100 5250 50  0001 C CNN
	1    10100 5250
	-1   0    0    -1  
$EndComp
Wire Wire Line
	9850 5250 9700 5250
Wire Wire Line
	9450 5250 9450 5350
Wire Wire Line
	9850 5350 9700 5350
Wire Wire Line
	8550 5550 9400 5550
Connection ~ 9700 5150
Wire Wire Line
	9700 5150 9850 5150
Connection ~ 9700 5250
Wire Wire Line
	9700 5250 9500 5250
Connection ~ 9700 5350
Wire Wire Line
	9700 5350 9600 5350
Wire Wire Line
	8300 5550 8550 5550
Wire Wire Line
	8400 5350 8850 5350
Wire Wire Line
	8250 5550 8550 5550
Wire Wire Line
	8250 5350 8850 5350
$Comp
L power:+5V #PWR029
U 1 1 6243E929
P 9700 4350
F 0 "#PWR029" H 9700 4200 50  0001 C CNN
F 1 "+5V" H 9715 4523 50  0000 C CNN
F 2 "" H 9700 4350 50  0001 C CNN
F 3 "" H 9700 4350 50  0001 C CNN
	1    9700 4350
	1    0    0    -1  
$EndComp
Wire Wire Line
	9700 4350 9700 4450
$Comp
L Device:C C25
U 1 1 6244180B
P 9150 4600
F 0 "C25" H 9265 4646 50  0000 L CNN
F 1 "100n" H 9265 4555 50  0000 L CNN
F 2 "Capacitor_THT:C_Rect_L7.0mm_W2.0mm_P5.00mm" H 9188 4450 50  0001 C CNN
F 3 "~" H 9150 4600 50  0001 C CNN
	1    9150 4600
	1    0    0    -1  
$EndComp
$Comp
L power:GNDREF #PWR027
U 1 1 6244771C
P 9150 4750
F 0 "#PWR027" H 9150 4500 50  0001 C CNN
F 1 "GNDREF" H 9155 4577 50  0001 C CNN
F 2 "" H 9150 4750 50  0001 C CNN
F 3 "" H 9150 4750 50  0001 C CNN
	1    9150 4750
	1    0    0    -1  
$EndComp
Wire Wire Line
	9150 4450 9600 4450
Connection ~ 9700 4450
Wire Wire Line
	9700 4450 9700 5050
Wire Wire Line
	8500 2600 8750 2600
Wire Wire Line
	9550 2550 9550 2600
Wire Wire Line
	9300 2600 9550 2600
Connection ~ 9550 2600
Wire Wire Line
	9550 2600 9550 2700
Wire Wire Line
	9550 2250 9900 2250
Wire Wire Line
	9550 2700 9900 2700
Wire Wire Line
	9900 2700 9900 2800
Connection ~ 9550 2700
Wire Wire Line
	9550 2700 9550 2750
$Comp
L power:GNDREF #PWR028
U 1 1 6245BA33
P 9550 3050
F 0 "#PWR028" H 9550 2800 50  0001 C CNN
F 1 "GNDREF" H 9555 2877 50  0001 C CNN
F 2 "" H 9550 3050 50  0001 C CNN
F 3 "" H 9550 3050 50  0001 C CNN
	1    9550 3050
	1    0    0    -1  
$EndComp
$Comp
L power:GNDREF #PWR030
U 1 1 6245BEC6
P 9900 2550
F 0 "#PWR030" H 9900 2300 50  0001 C CNN
F 1 "GNDREF" H 9905 2377 50  0001 C CNN
F 2 "" H 9900 2550 50  0001 C CNN
F 3 "" H 9900 2550 50  0001 C CNN
	1    9900 2550
	1    0    0    -1  
$EndComp
$Comp
L power:GNDREF #PWR031
U 1 1 6245C334
P 9900 3100
F 0 "#PWR031" H 9900 2850 50  0001 C CNN
F 1 "GNDREF" H 9905 2927 50  0001 C CNN
F 2 "" H 9900 3100 50  0001 C CNN
F 3 "" H 9900 3100 50  0001 C CNN
	1    9900 3100
	1    0    0    -1  
$EndComp
Wire Wire Line
	9650 1850 8750 1850
Wire Wire Line
	8750 1850 8750 2600
Connection ~ 8750 2600
Wire Wire Line
	8750 2600 9000 2600
Text GLabel 10500 1850 0    50   Input ~ 0
TO_RECEIVER
$Comp
L Device:C C2
U 1 1 624874F7
P 850 2950
F 0 "C2" H 965 2996 50  0000 L CNN
F 1 "100n" H 965 2905 50  0000 L CNN
F 2 "Capacitor_THT:C_Rect_L7.0mm_W2.0mm_P5.00mm" H 888 2800 50  0001 C CNN
F 3 "~" H 850 2950 50  0001 C CNN
	1    850  2950
	1    0    0    1   
$EndComp
$Comp
L Device:C C3
U 1 1 624874FD
P 1550 1650
F 0 "C3" H 1665 1696 50  0000 L CNN
F 1 "100n" H 1665 1605 50  0000 L CNN
F 2 "Capacitor_THT:C_Rect_L7.0mm_W2.0mm_P5.00mm" H 1588 1500 50  0001 C CNN
F 3 "~" H 1550 1650 50  0001 C CNN
	1    1550 1650
	0    -1   -1   0   
$EndComp
$Comp
L Device:C C1
U 1 1 62487503
P 850 2400
F 0 "C1" H 965 2446 50  0000 L CNN
F 1 "100n" H 965 2355 50  0000 L CNN
F 2 "Capacitor_THT:C_Rect_L7.0mm_W2.0mm_P5.00mm" H 888 2250 50  0001 C CNN
F 3 "~" H 850 2400 50  0001 C CNN
	1    850  2400
	1    0    0    1   
$EndComp
$Comp
L Device:R R2
U 1 1 62487509
P 1200 2900
F 0 "R2" H 1270 2946 50  0000 L CNN
F 1 "15k" H 1270 2855 50  0000 L CNN
F 2 "Resistor_THT:R_Axial_DIN0207_L6.3mm_D2.5mm_P2.54mm_Vertical" V 1130 2900 50  0001 C CNN
F 3 "~" H 1200 2900 50  0001 C CNN
	1    1200 2900
	1    0    0    1   
$EndComp
$Comp
L Device:R R1
U 1 1 6248750F
P 1200 2400
F 0 "R1" H 1270 2446 50  0000 L CNN
F 1 "15k" H 1270 2355 50  0000 L CNN
F 2 "Resistor_THT:R_Axial_DIN0207_L6.3mm_D2.5mm_P2.54mm_Vertical" V 1130 2400 50  0001 C CNN
F 3 "~" H 1200 2400 50  0001 C CNN
	1    1200 2400
	1    0    0    1   
$EndComp
$Comp
L Device:R R3
U 1 1 62487515
P 1600 2600
F 0 "R3" H 1670 2646 50  0000 L CNN
F 1 "10k" H 1670 2555 50  0000 L CNN
F 2 "Resistor_THT:R_Axial_DIN0207_L6.3mm_D2.5mm_P2.54mm_Vertical" V 1530 2600 50  0001 C CNN
F 3 "~" H 1600 2600 50  0001 C CNN
	1    1600 2600
	0    -1   1    0   
$EndComp
Text GLabel 2150 3050 0    50   Input ~ 0
S2
Text GLabel 2150 2900 0    50   Input ~ 0
S1
Text GLabel 2150 2750 0    50   Input ~ 0
S0
Wire Wire Line
	2300 2800 2150 2800
Wire Wire Line
	2150 2800 2150 2750
Wire Wire Line
	2300 2900 2150 2900
Wire Wire Line
	2300 3000 2150 3000
Wire Wire Line
	2150 3000 2150 3050
Wire Wire Line
	1200 2550 1200 2600
Wire Wire Line
	1450 2600 1200 2600
Connection ~ 1200 2600
Wire Wire Line
	1200 2600 1200 2700
Wire Wire Line
	1200 2250 850  2250
Wire Wire Line
	1200 2700 850  2700
Wire Wire Line
	850  2700 850  2800
Connection ~ 1200 2700
Wire Wire Line
	1200 2700 1200 2750
$Comp
L power:GNDREF #PWR03
U 1 1 6248752D
P 1200 3050
F 0 "#PWR03" H 1200 2800 50  0001 C CNN
F 1 "GNDREF" H 1205 2877 50  0001 C CNN
F 2 "" H 1200 3050 50  0001 C CNN
F 3 "" H 1200 3050 50  0001 C CNN
	1    1200 3050
	-1   0    0    -1  
$EndComp
$Comp
L power:GNDREF #PWR01
U 1 1 62487533
P 850 2550
F 0 "#PWR01" H 850 2300 50  0001 C CNN
F 1 "GNDREF" H 855 2377 50  0001 C CNN
F 2 "" H 850 2550 50  0001 C CNN
F 3 "" H 850 2550 50  0001 C CNN
	1    850  2550
	-1   0    0    -1  
$EndComp
$Comp
L power:GNDREF #PWR02
U 1 1 62487539
P 850 3100
F 0 "#PWR02" H 850 2850 50  0001 C CNN
F 1 "GNDREF" H 855 2927 50  0001 C CNN
F 2 "" H 850 3100 50  0001 C CNN
F 3 "" H 850 3100 50  0001 C CNN
	1    850  3100
	-1   0    0    -1  
$EndComp
Wire Wire Line
	1750 2600 2000 2600
Wire Wire Line
	1700 1650 2000 1650
Wire Wire Line
	2000 1650 2000 2600
Connection ~ 2000 2600
Wire Wire Line
	2000 2600 2300 2600
$Comp
L Device:R R5
U 1 1 624EDC6D
P 8700 3450
F 0 "R5" H 8770 3496 50  0000 L CNN
F 1 "1k" H 8770 3405 50  0000 L CNN
F 2 "Resistor_THT:R_Axial_DIN0207_L6.3mm_D2.5mm_P2.54mm_Vertical" V 8630 3450 50  0001 C CNN
F 3 "~" H 8700 3450 50  0001 C CNN
	1    8700 3450
	1    0    0    -1  
$EndComp
Wire Wire Line
	8700 3300 8700 3200
Wire Wire Line
	8700 3200 8500 3200
$Comp
L power:GNDREF #PWR025
U 1 1 624FF1FE
P 8700 3600
F 0 "#PWR025" H 8700 3350 50  0001 C CNN
F 1 "GNDREF" H 8705 3427 50  0001 C CNN
F 2 "" H 8700 3600 50  0001 C CNN
F 3 "" H 8700 3600 50  0001 C CNN
	1    8700 3600
	1    0    0    -1  
$EndComp
Text GLabel 700  1450 2    50   Output ~ 0
FROM_ANTENNA
$Comp
L Device:R R4
U 1 1 62550A52
P 2000 3450
F 0 "R4" H 2070 3496 50  0000 L CNN
F 1 "1k" H 2070 3405 50  0000 L CNN
F 2 "Resistor_THT:R_Axial_DIN0207_L6.3mm_D2.5mm_P2.54mm_Vertical" V 1930 3450 50  0001 C CNN
F 3 "~" H 2000 3450 50  0001 C CNN
	1    2000 3450
	1    0    0    -1  
$EndComp
Wire Wire Line
	2300 3200 2000 3200
Wire Wire Line
	2000 3200 2000 3300
$Comp
L power:GNDREF #PWR04
U 1 1 62555FAA
P 2000 3600
F 0 "#PWR04" H 2000 3350 50  0001 C CNN
F 1 "GNDREF" H 2005 3427 50  0001 C CNN
F 2 "" H 2000 3600 50  0001 C CNN
F 3 "" H 2000 3600 50  0001 C CNN
	1    2000 3600
	1    0    0    -1  
$EndComp
$Comp
L power:GNDREF #PWR05
U 1 1 6258488C
P 2600 3500
F 0 "#PWR05" H 2600 3250 50  0001 C CNN
F 1 "GNDREF" H 2605 3327 50  0001 C CNN
F 2 "" H 2600 3500 50  0001 C CNN
F 3 "" H 2600 3500 50  0001 C CNN
	1    2600 3500
	1    0    0    -1  
$EndComp
$Comp
L power:GNDREF #PWR06
U 1 1 62585DF3
P 2700 3500
F 0 "#PWR06" H 2700 3250 50  0001 C CNN
F 1 "GNDREF" H 2705 3327 50  0001 C CNN
F 2 "" H 2700 3500 50  0001 C CNN
F 3 "" H 2700 3500 50  0001 C CNN
	1    2700 3500
	1    0    0    -1  
$EndComp
Text Label 7700 7500 0    50   ~ 0
HF-Filter
Text Label 8200 7500 0    50   ~ 0
IK8YFW-2022
Wire Wire Line
	1400 1650 1400 1450
Wire Wire Line
	1400 1450 700  1450
$Comp
L power:+5V #PWR0101
U 1 1 623832AE
P 8150 1850
F 0 "#PWR0101" H 8150 1700 50  0001 C CNN
F 1 "+5V" H 8165 2023 50  0000 C CNN
F 2 "" H 8150 1850 50  0001 C CNN
F 3 "" H 8150 1850 50  0001 C CNN
	1    8150 1850
	1    0    0    -1  
$EndComp
Wire Wire Line
	8200 2400 8200 1850
Wire Wire Line
	8200 1850 8150 1850
Connection ~ 8150 1850
Wire Wire Line
	8150 1850 8100 1850
$Comp
L power:+5V #PWR0102
U 1 1 6238DBD0
P 2550 2050
F 0 "#PWR0102" H 2550 1900 50  0001 C CNN
F 1 "+5V" H 2565 2223 50  0000 C CNN
F 2 "" H 2550 2050 50  0001 C CNN
F 3 "" H 2550 2050 50  0001 C CNN
	1    2550 2050
	1    0    0    -1  
$EndComp
Wire Wire Line
	2600 2400 2600 2050
Wire Wire Line
	2600 2050 2550 2050
$Comp
L Connector:Screw_Terminal_01x02 J1
U 1 1 624CF3BD
P 1050 3750
F 0 "J1" H 1130 3742 50  0000 L CNN
F 1 "Screw_Terminal" H 1130 3651 50  0000 L CNN
F 2 "" H 1050 3750 50  0001 C CNN
F 3 "~" H 1050 3750 50  0001 C CNN
	1    1050 3750
	1    0    0    -1  
$EndComp
$Comp
L Connector:Screw_Terminal_01x04 J2
U 1 1 624D1792
P 9900 6000
F 0 "J2" H 9980 5992 50  0000 L CNN
F 1 "Screw_Terminal" H 9980 5901 50  0000 L CNN
F 2 "" H 9900 6000 50  0001 C CNN
F 3 "~" H 9900 6000 50  0001 C CNN
	1    9900 6000
	1    0    0    -1  
$EndComp
Wire Wire Line
	9700 5900 9700 5450
Wire Wire Line
	9700 5250 9700 5150
Wire Wire Line
	9700 5350 9700 5250
Connection ~ 9700 5450
Wire Wire Line
	9700 5450 9700 5350
Wire Wire Line
	9700 6000 9500 6000
Wire Wire Line
	9500 6000 9500 5250
Connection ~ 9500 5250
Wire Wire Line
	9500 5250 9450 5250
Wire Wire Line
	9700 6100 9400 6100
Wire Wire Line
	9400 6100 9400 5550
Wire Wire Line
	9350 5150 9700 5150
Wire Wire Line
	9150 5150 9700 5150
Connection ~ 9400 5550
Wire Wire Line
	9400 5550 9600 5550
Wire Wire Line
	9700 6200 9600 6200
Wire Wire Line
	9600 4450 9600 5350
Connection ~ 9600 4450
Wire Wire Line
	9600 4450 9700 4450
Connection ~ 9600 5350
Wire Wire Line
	9600 5350 9600 5550
Connection ~ 9600 5550
Wire Wire Line
	9600 5550 9600 6200
$Comp
L power:GNDREF #PWR0103
U 1 1 624E9ED6
P 850 3950
F 0 "#PWR0103" H 850 3700 50  0001 C CNN
F 1 "GNDREF" H 855 3777 50  0001 C CNN
F 2 "" H 850 3950 50  0001 C CNN
F 3 "" H 850 3950 50  0001 C CNN
	1    850  3950
	1    0    0    -1  
$EndComp
$Comp
L power:+5V #PWR0104
U 1 1 624EAE5B
P 750 3650
F 0 "#PWR0104" H 750 3500 50  0001 C CNN
F 1 "+5V" H 765 3823 50  0000 C CNN
F 2 "" H 750 3650 50  0001 C CNN
F 3 "" H 750 3650 50  0001 C CNN
	1    750  3650
	1    0    0    -1  
$EndComp
Wire Wire Line
	850  3750 750  3750
Wire Wire Line
	750  3750 750  3650
Wire Wire Line
	850  3850 850  3950
Connection ~ 850  3950
Wire Wire Line
	850  3950 850  4000
$EndSCHEMATC