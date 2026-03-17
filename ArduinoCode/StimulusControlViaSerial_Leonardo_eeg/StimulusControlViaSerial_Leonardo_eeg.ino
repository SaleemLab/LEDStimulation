#include <math.h>  //

#define CLOCK_FREQ 16000000                // Arduino Leonardo clock frequency (16 MHz)
#define TABLE_SIZE 256                     // Number of samples in the wavetable
#define FWN_TABLE_SIZE 375                 // Number of samples in the wavetable
#define TARGET_RECONFIG_INTERVAL_US 1000L  // FOR FREQUENCY SWEEP 1000 microseconds = 1 millisecond

volatile unsigned long printSequenceNum = 0;  // for serial debugging

// --- Faster PRNG (xorshift32) ---
uint32_t xorshift32_state = 1;  // Seed with a non-zero value.

uint32_t fast_rand32() {
  xorshift32_state ^= (xorshift32_state << 13);
  xorshift32_state ^= (xorshift32_state >> 17);
  xorshift32_state ^= (xorshift32_state << 5);
  return xorshift32_state;
}



// gamma-correction LUTS
const uint16_t PROGMEM ChA1LUT[1041] = {
  0, 6, 11, 14, 16, 19, 21, 23, 25, 27,
  28, 30, 32, 34, 35, 37, 38, 40, 41, 43,
  44, 45, 46, 47, 49, 50, 51, 52, 54, 55,
  56, 57, 59, 60, 61, 62, 63, 64, 66, 67,
  68, 69, 70, 71, 72, 73, 74, 76, 77, 78,
  79, 80, 82, 83, 84, 85, 86, 87, 89, 90,
  91, 92, 93, 94, 95, 96, 97, 99, 100, 101,
  102, 103, 104, 105, 106, 107, 108, 109, 111, 112,
  113, 114, 115, 116, 117, 118, 119, 120, 121, 122,
  123, 124, 125, 126, 128, 129, 130, 131, 132, 133,
  134, 135, 137, 138, 139, 140, 141, 142, 143, 145,
  146, 147, 148, 149, 150, 151, 152, 153, 154, 155,
  156, 157, 158, 159, 160, 161, 162, 164, 165, 166,
  167, 168, 169, 170, 171, 172, 173, 174, 175, 176,
  177, 178, 179, 180, 182, 183, 184, 185, 186, 187,
  188, 189, 190, 191, 192, 193, 194, 195, 196, 197,
  198, 199, 200, 202, 203, 204, 205, 206, 207, 208,
  209, 210, 211, 213, 214, 215, 216, 217, 218, 219,
  220, 221, 222, 223, 224, 225, 226, 227, 228, 229,
  230, 231, 232, 233, 235, 236, 237, 238, 239, 240,
  241, 242, 243, 244, 245, 246, 247, 248, 249, 250,
  251, 252, 253, 254, 255, 257, 258, 259, 260, 261,
  262, 263, 264, 265, 266, 267, 268, 269, 271, 272,
  273, 274, 275, 276, 277, 278, 279, 280, 281, 282,
  283, 284, 285, 286, 287, 288, 289, 290, 291, 292,
  294, 295, 296, 297, 298, 299, 300, 301, 302, 303,
  304, 305, 306, 307, 308, 309, 310, 311, 312, 314,
  315, 316, 317, 318, 319, 320, 321, 322, 323, 324,
  325, 326, 327, 328, 329, 330, 331, 332, 333, 334,
  335, 336, 337, 338, 339, 340, 341, 342, 343, 345,
  346, 347, 348, 349, 350, 351, 352, 353, 354, 355,
  356, 357, 358, 359, 360, 361, 362, 363, 364, 364,
  366, 367, 368, 369, 370, 371, 372, 373, 374, 375,
  376, 377, 378, 379, 380, 381, 382, 383, 384, 385,
  386, 387, 388, 389, 390, 391, 392, 393, 394, 395,
  396, 397, 398, 399, 400, 401, 402, 403, 404, 405,
  406, 407, 408, 409, 410, 411, 412, 413, 414, 415,
  416, 417, 418, 419, 420, 421, 422, 424, 425, 426,
  427, 428, 428, 429, 430, 431, 432, 433, 434, 435,
  436, 437, 438, 439, 440, 441, 442, 443, 444, 445,
  446, 447, 448, 449, 450, 451, 452, 453, 454, 455,
  456, 457, 458, 459, 460, 461, 462, 463, 464, 465,
  466, 466, 467, 468, 469, 470, 471, 472, 473, 474,
  476, 477, 478, 479, 480, 480, 481, 482, 483, 484,
  485, 486, 487, 488, 489, 490, 491, 492, 493, 494,
  495, 496, 497, 498, 499, 500, 501, 502, 503, 504,
  505, 506, 507, 508, 509, 510, 511, 512, 513, 513,
  514, 515, 516, 517, 518, 519, 520, 521, 522, 523,
  524, 525, 526, 527, 528, 529, 530, 531, 532, 533,
  534, 534, 535, 536, 537, 538, 539, 540, 541, 542,
  543, 544, 545, 546, 547, 548, 549, 550, 551, 552,
  553, 554, 555, 556, 557, 558, 559, 560, 561, 562,
  563, 564, 565, 566, 566, 567, 568, 569, 570, 571,
  572, 573, 574, 575, 576, 577, 578, 579, 580, 581,
  582, 583, 584, 585, 585, 586, 587, 588, 589, 590,
  591, 592, 593, 594, 595, 595, 596, 597, 598, 599,
  600, 601, 602, 603, 604, 605, 606, 607, 608, 609,
  610, 611, 612, 613, 614, 614, 615, 616, 617, 618,
  619, 620, 621, 622, 623, 624, 625, 626, 627, 628,
  629, 630, 631, 632, 633, 634, 635, 635, 636, 637,
  638, 639, 640, 641, 642, 642, 643, 644, 645, 646,
  647, 648, 649, 650, 651, 652, 653, 654, 655, 656,
  657, 658, 659, 660, 661, 662, 663, 664, 665, 666,
  667, 667, 668, 669, 670, 671, 672, 673, 673, 674,
  675, 676, 677, 678, 679, 680, 681, 682, 683, 684,
  685, 686, 686, 687, 688, 689, 690, 691, 692, 693,
  694, 694, 695, 696, 697, 698, 699, 700, 701, 702,
  703, 704, 705, 706, 707, 707, 708, 709, 710, 711,
  712, 713, 714, 715, 716, 717, 718, 719, 720, 721,
  722, 723, 724, 724, 725, 726, 727, 728, 729, 730,
  731, 732, 733, 734, 735, 736, 737, 738, 739, 740,
  741, 742, 743, 743, 744, 745, 746, 747, 748, 749,
  749, 750, 751, 752, 753, 754, 755, 756, 757, 758,
  759, 760, 761, 762, 763, 764, 765, 766, 767, 768,
  769, 770, 770, 771, 772, 773, 774, 775, 775, 776,
  777, 778, 779, 779, 780, 781, 782, 783, 784, 785,
  786, 787, 788, 789, 790, 791, 792, 793, 794, 795,
  795, 796, 797, 798, 799, 800, 801, 802, 803, 804,
  805, 806, 806, 807, 808, 809, 810, 811, 812, 813,
  814, 815, 816, 817, 817, 818, 819, 820, 821, 822,
  823, 824, 824, 825, 826, 827, 828, 829, 830, 830,
  831, 832, 833, 834, 835, 836, 837, 838, 839, 840,
  840, 841, 842, 843, 844, 845, 846, 846, 847, 848,
  849, 850, 851, 851, 852, 853, 854, 855, 856, 857,
  858, 859, 860, 861, 862, 863, 863, 864, 865, 866,
  867, 868, 869, 870, 871, 872, 873, 873, 874, 875,
  876, 877, 877, 878, 879, 880, 880, 881, 882, 883,
  883, 884, 885, 886, 887, 888, 889, 890, 891, 892,
  893, 893, 894, 895, 896, 897, 898, 899, 899, 900,
  901, 902, 903, 904, 904, 905, 906, 907, 908, 909,
  910, 911, 912, 912, 913, 914, 915, 916, 917, 918,
  919, 919, 920, 921, 922, 923, 924, 925, 925, 926,
  927, 928, 929, 930, 931, 931, 932, 933, 934, 935,
  936, 937, 938, 938, 939, 940, 941, 942, 943, 944,
  945, 946, 946, 947, 948, 949, 950, 950, 951, 952,
  953, 953, 954, 955, 956, 956, 957, 958, 959, 960,
  961, 962, 963, 964, 964, 965, 966, 967, 968, 969,
  970, 971, 972, 973, 974, 974, 975, 976, 977, 978,
  979, 980, 980, 981, 982, 983, 983, 984, 985, 986,
  987, 987, 988, 989, 990, 991, 992, 993, 994, 994,
  995, 996, 997, 998, 999, 1000, 1000, 1001, 1002, 1003,
  1004, 1004, 1005, 1006, 1007, 1007, 1008, 1009, 1010, 1011,
  1012, 1013, 1014, 1014, 1015, 1016, 1017, 1018, 1019, 1020,
  1021, 1022, 1022, 1023, 1024, 1025, 1026, 1027, 1028, 1029,
  1030
};

const uint16_t PROGMEM ChA2LUT[1041] = {
  0, 45, 49, 52, 54, 56, 59, 61, 63, 65,
  68, 70, 73, 74, 76, 77, 79, 80, 82, 83,
  85, 87, 90, 92, 94, 96, 98, 100, 102, 104,
  106, 107, 109, 110, 112, 113, 114, 116, 117, 119,
  120, 122, 123, 125, 127, 128, 130, 131, 133, 134,
  136, 137, 139, 141, 142, 144, 145, 146, 147, 148,
  149, 150, 151, 153, 154, 155, 156, 157, 159, 161,
  162, 164, 166, 167, 168, 169, 170, 171, 172, 173,
  175, 176, 177, 178, 179, 180, 181, 182, 183, 184,
  185, 186, 188, 189, 190, 192, 193, 194, 196, 197,
  198, 199, 200, 201, 202, 203, 204, 205, 206, 207,
  208, 210, 211, 212, 213, 214, 215, 216, 217, 218,
  219, 220, 222, 223, 224, 225, 226, 227, 228, 229,
  230, 231, 233, 234, 235, 236, 237, 238, 239, 240,
  241, 242, 243, 244, 245, 246, 247, 249, 250, 250,
  251, 252, 253, 254, 255, 255, 256, 257, 258, 259,
  259, 260, 261, 263, 264, 265, 266, 267, 268, 269,
  270, 271, 272, 273, 274, 275, 276, 277, 278, 279,
  280, 281, 282, 283, 284, 285, 286, 287, 288, 289,
  290, 291, 292, 294, 295, 296, 297, 298, 299, 300,
  301, 302, 303, 304, 305, 306, 307, 307, 308, 309,
  310, 311, 312, 313, 314, 315, 316, 316, 317, 318,
  319, 320, 321, 322, 323, 324, 325, 326, 327, 328,
  329, 330, 331, 332, 333, 334, 335, 336, 337, 339,
  340, 341, 342, 343, 344, 345, 345, 346, 347, 347,
  348, 349, 350, 350, 351, 352, 352, 353, 354, 355,
  356, 357, 359, 360, 361, 362, 363, 364, 365, 366,
  367, 368, 369, 370, 371, 373, 374, 375, 375, 376,
  377, 377, 378, 379, 380, 380, 381, 382, 382, 383,
  384, 384, 385, 386, 387, 388, 389, 390, 391, 392,
  393, 394, 395, 396, 396, 397, 398, 399, 400, 401,
  402, 403, 404, 405, 406, 407, 408, 409, 410, 411,
  412, 413, 414, 415, 416, 417, 418, 419, 420, 421,
  422, 423, 424, 425, 426, 427, 428, 428, 429, 430,
  430, 431, 432, 432, 433, 434, 434, 435, 436, 436,
  437, 438, 440, 441, 442, 443, 444, 445, 446, 447,
  448, 449, 450, 451, 452, 453, 454, 455, 455, 456,
  457, 458, 459, 460, 461, 462, 463, 464, 465, 466,
  467, 468, 469, 469, 470, 471, 472, 473, 474, 475,
  476, 477, 478, 479, 479, 480, 481, 482, 482, 483,
  484, 485, 486, 486, 487, 488, 489, 490, 490, 491,
  492, 493, 494, 495, 496, 497, 498, 498, 499, 500,
  501, 502, 503, 504, 504, 505, 506, 507, 508, 509,
  510, 511, 512, 513, 514, 515, 516, 517, 518, 519,
  520, 521, 521, 522, 523, 523, 524, 525, 526, 526,
  527, 528, 529, 529, 530, 531, 532, 532, 533, 534,
  535, 536, 537, 537, 538, 539, 540, 541, 542, 543,
  544, 545, 546, 546, 547, 548, 549, 550, 551, 552,
  553, 554, 555, 556, 557, 558, 558, 559, 560, 561,
  562, 563, 564, 565, 565, 566, 567, 568, 569, 570,
  571, 572, 573, 574, 575, 576, 577, 578, 579, 580,
  581, 582, 583, 584, 584, 585, 586, 586, 587, 588,
  588, 589, 590, 590, 591, 592, 592, 593, 594, 595,
  596, 597, 597, 598, 599, 600, 601, 602, 603, 604,
  605, 606, 606, 607, 608, 609, 610, 611, 612, 613,
  614, 614, 615, 616, 617, 618, 619, 620, 621, 622,
  623, 624, 625, 626, 627, 627, 628, 629, 630, 631,
  631, 632, 633, 634, 634, 635, 636, 637, 637, 638,
  639, 640, 640, 641, 642, 643, 643, 644, 645, 646,
  647, 648, 649, 650, 651, 652, 653, 654, 655, 656,
  657, 658, 659, 660, 661, 662, 663, 664, 665, 666,
  666, 667, 668, 668, 669, 670, 670, 671, 672, 672,
  673, 674, 674, 675, 676, 676, 677, 678, 679, 680,
  681, 683, 684, 685, 686, 687, 687, 688, 689, 690,
  690, 691, 692, 693, 694, 694, 695, 696, 697, 697,
  698, 699, 700, 701, 702, 703, 703, 704, 705, 706,
  707, 708, 709, 710, 711, 711, 712, 713, 714, 715,
  716, 717, 718, 719, 720, 720, 721, 722, 723, 724,
  725, 726, 726, 727, 728, 729, 730, 731, 731, 732,
  733, 734, 735, 736, 736, 737, 738, 739, 740, 741,
  742, 742, 743, 744, 745, 746, 747, 748, 748, 749,
  750, 751, 752, 752, 753, 754, 755, 756, 757, 757,
  758, 759, 760, 761, 762, 763, 764, 765, 766, 767,
  768, 769, 770, 771, 771, 772, 773, 774, 774, 775,
  776, 776, 777, 778, 778, 779, 780, 781, 782, 783,
  784, 785, 786, 787, 788, 789, 790, 791, 791, 792,
  793, 794, 795, 795, 796, 797, 798, 798, 799, 800,
  801, 801, 802, 803, 804, 805, 805, 806, 807, 808,
  809, 809, 810, 811, 812, 813, 814, 814, 815, 816,
  817, 818, 819, 820, 821, 822, 823, 824, 825, 826,
  827, 828, 829, 830, 831, 832, 833, 834, 834, 835,
  836, 836, 837, 838, 838, 839, 840, 840, 841, 842,
  842, 844, 845, 846, 847, 849, 850, 851, 852, 853,
  854, 854, 855, 856, 856, 857, 857, 858, 859, 859,
  860, 860, 861, 862, 862, 863, 864, 866, 867, 868,
  870, 871, 873, 874, 875, 876, 876, 877, 878, 879,
  880, 880, 881, 882, 883, 884, 884, 885, 886, 886,
  887, 888, 888, 889, 890, 891, 891, 892, 893, 893,
  894, 895, 896, 896, 897, 898, 899, 900, 900, 901,
  902, 903, 904, 904, 905, 906, 907, 908, 909, 910,
  910, 911, 912, 913, 914, 915, 916, 917, 918, 919,
  920, 921, 922, 923, 924, 925, 926, 927, 928, 928,
  929, 930, 930, 931, 932, 933, 933, 934, 935, 935,
  936, 937, 938, 939, 940, 940, 941, 942, 943, 944,
  945, 945, 946, 947, 948, 948, 949, 950, 950, 951,
  952, 952, 953, 954, 954, 955, 956, 956, 957, 958,
  959, 960, 961, 962, 963, 964, 965, 966, 967, 968,
  969, 971, 972, 973, 974, 975, 976, 977, 978, 979,
  979, 980, 980, 981, 982, 982, 983, 983, 984, 985,
  985, 986, 986, 987, 988, 989, 990, 991, 992, 993,
  994, 995, 996, 997, 998, 999, 1000, 1000, 1001, 1002,
  1003, 1003, 1004, 1005, 1006, 1007, 1007, 1008, 1009, 1010,
  1011, 1011, 1012, 1013, 1014, 1014, 1015, 1016, 1017, 1017,
  1018, 1019, 1020, 1021, 1022, 1024, 1025, 1026, 1027, 1028,
  1030
};


const uint16_t *currentChALUT;
const uint16_t *currentChBLUT;


// PWM variables
long prescaler;
uint16_t TOP;  // set by the clock
long TopLumi;  // use to limit max luminance
long MidLumi;
long desiredPWMFrequency = 7680;  // 2000?
//float dutyCycle;

// channel selection
bool useChA = true;
bool useChB = true;

// serial
const byte numChars = 50;
char receivedChars[numChars];  // an array to store the received data
bool newData = false;

// stimulus selection char
String FirstChar;

// Array to hold the wavetable
uint16_t sineWaveTable[TABLE_SIZE];

// wavetable step size
volatile int stepSize = 1;
volatile int tableIndexA = 0;  //for sinewave table
volatile int tableIndexB = 0;  //for sinewave table
volatile float contrastMultA;
volatile float contrastMultB;
volatile int tableEnvIndex = 0;  // for sinewave table contrast envelope
volatile int envCount = 0;       // contrast envelope counter
// contrast-envelope counter and contrast multiplier
volatile int nEnvCounts = 1;
volatile float contrastMult = 0;
volatile uint16_t contrastMultInt = 0;
volatile unsigned int completedCycles = 0; // Tracks full sine wave cycles

// white noise parameters
volatile float target_mean;  // Desired mean of the final output
volatile float target_std;   // Desired standard deviation of the final output

volatile uint16_t finalRandNumber_A = 0;  // Final output value 1
volatile uint16_t finalRandNumber_B = 0;  // Final output value 2
volatile float map_float_min;
volatile float map_float_max;
volatile float mapped_float_A;
volatile float mapped_float_B;

const int CLT_N = 10;                 // Number of uniform samples for Central Limit Theorem
                                      // Higher N => better Gaussian approximation, but slower.
                                      // (N >= 10 or 12 is common)
const int RANDOM_UPPER_BOUND = 1024;  // The argument 'M' for random(M). Generates [0, M-1]. Now using fast_rand32(), so must be power of 2 (e.g. 1024)
// --- Derived Constants (calculated once for efficiency) ---
// Pre-calculate values needed for N(0,1) generation based on configuration
const float UNIFORM_MEAN = (float)(RANDOM_UPPER_BOUND - 1) / 2.0;
const float SUM_EXPECTED_MEAN = (float)CLT_N * UNIFORM_MEAN;
// Variance of U[0, M-1] = (M^2 - 1) / 12
const double UNIFORM_VARIANCE = (pow((double)RANDOM_UPPER_BOUND, 2) - 1.0) / 12.0;
// Variance of Sum = N * Var(U)
const double SUM_VARIANCE = (double)CLT_N * UNIFORM_VARIANCE;
// Std Dev of Sum = sqrt(Var(Sum))
const double SUM_STD_DEV = sqrt(SUM_VARIANCE);
// Scale factor to achieve N(0,1) = 1 / Std Dev of Sum
const float NORMALIZE_SCALE_FACTOR = (SUM_STD_DEV > 0) ? (1.0 / (float)SUM_STD_DEV) : 0.0;

//update time and PWM value
volatile long updateTime;
volatile long randNumber = 0;
volatile int currentDist;

// array for frozen white noise values
volatile uint16_t frozenWhiteNoiseTable[FWN_TABLE_SIZE];  // max number of values (8ms update time for a 3s stimulus)
volatile int tableIndexFWN = 0;                           //for sinewave table

// flicker state
volatile bool toggleState = false;  // Flag to track the current state


void setup() {
  pinMode(9, OUTPUT);   // Pin 9 controlled by Timer1 (Channel A)
  pinMode(10, OUTPUT);  // Pin 9 controlled by Timer1 (Channel B)

  pinMode(4, OUTPUT);  // Pin 4 indicator pin, e.g. for sinewave cycles, changes in duty cycle during white noise, changes in contrast etc.
  pinMode(5, OUTPUT);  // Pin 5 stim ON or OFF pin

  PORTD &= ~(1 << PIND4);   // Ensure Pin 4 is set to LOW by changing register directly
  PORTC &= ~(1 << PORTC6);  // Ensure Pin 5 is set to LOW

  // Check fixed derived constants for white noise generation
  if (NORMALIZE_SCALE_FACTOR == 0.0) {
    Serial.println("ERROR: Normalization scale factor is zero. Halting.");
    while (1)
      ;
  }

  // set default LUTs to "1"
  currentChALUT = ChA1LUT;
  currentChBLUT = ChA1LUT;


  Serial.begin(115200);


  // Constrain the desired PWM frequency to be a multiple of TABLE_SIZE
  // (probably not important since precision seems low at high frequenies)
  desiredPWMFrequency = constrainFrequency(desiredPWMFrequency);

  // Calculate the prescaler and TOP value for the constrained frequency
  TOP = calculatePrescalerAndTOP(desiredPWMFrequency, prescaler);

  // Apply the prescaler and TOP value to Timer1
  configureTimer1(prescaler, TOP);

  TopLumi = TOP;  // set default TopLumi (max duty cycle) as the actual TOP (i.e. 100% for now)
  MidLumi = TOP / 2;

  // Generate the sine wave LUT based on the Timer1 config and TopLumi
  generateSineWaveTable(TopLumi);

  // initialise random number to 50% duty cyle for white noise stimuli
  randomSeed(0);
  finalRandNumber_A = TopLumi / 2;
  finalRandNumber_B = TopLumi / 2;



  if (useChA) { setChA(TopLumi / 2); }  // Set pin 9 to 50% duty cycle as default
  if (useChB) { setChB(TopLumi / 2); }  // Set pin 10 to 50% duty cycle as default

  //delay(5000);
  //whiteNoise(10000, 10);
  //outputSinewave(10,10000);

  //ToggleLEDTest(8000);
}

// loop runs checking for new serial input
void loop() {
  GetSerialInput();
  if (newData) {
    newData = false;
    ActionSerial();
  }
}



//////////////////////////////////////// HANDLE SERIAL INPUT //////////////////////////////////////
void GetSerialInput() {  // part of code taken from http://forum.arduino.cc/index.php?topic=396450.0
  static byte ndx = 0;
  char endMarker = '\r';
  char rc;
  if (Serial.available() > 0) {
    rc = Serial.read();
    if (rc != endMarker) {
      receivedChars[ndx] = rc;
      ndx++;
      if (ndx >= numChars) {
        ndx = numChars - 1;
      }
    } else {                      // serial message finished
      receivedChars[ndx] = '\0';  // terminate the string
      ndx = 0;
      newData = true;
    }
  }
}

void ActionSerial() {  // Actions serial data by choosing appropriate stimulation
  Serial.print("rc: ");
  Serial.print(receivedChars);
  Serial.print("\n");
  char delimiters[] = ",";
  char *token;
  uint8_t idx = 0;
#define MAX_VALS 50  // max required? freq, duration, contrast, carrier freq?
  char *serialVals[MAX_VALS];
  token = strtok(receivedChars, ",");


  while (token != NULL) {
    //Serial.println( token );
    if (idx < MAX_VALS)
      serialVals[idx++] = token;
    token = strtok(NULL, ",");
  }

  FirstChar = serialVals[0];

  if (FirstChar == "s")  // sinusoidal flicker
  {
    long stimulusDuration = atof(serialVals[1]);
    float frequency = atof(serialVals[2]);
    float phaseA = atof(serialVals[3]);
    float phaseB = atof(serialVals[4]);
    float contrastA = atof(serialVals[5]);
    float contrastB = atof(serialVals[6]);

    //Serial.println("Stim: Sinusoidal dimming");
    //Serial.flush();
    //Serial.print("Stim duration: ");
    //Serial.println(stimulusDuration);
    //Serial.flush();
    //Serial.print("Frequency: ");
    //Serial.println(frequency);
    //Serial.flush();

    outputSinewave(frequency, stimulusDuration, phaseA, phaseB, contrastA, contrastB);

  } else if (FirstChar == "wn")  // white noise
  {
    long stimulusDuration = atof(serialVals[1]);
    long updateTime = atof(serialVals[2]);
    float frac_target_mean = atof(serialVals[3]);
    float frac_target_std = atof(serialVals[4]);
    //Serial.println("Stim: White noise");
    //Serial.print("Stim duration: ");
    // Serial.println(stimulusDuration);
    //Serial.flush();
    // Serial.print("Update time: ");
    // Serial.println(updateTime);
    //  Serial.flush();
    //Serial.println(stimulusDuration);
    whiteNoise(updateTime, stimulusDuration, frac_target_mean, frac_target_std);

  } else if (FirstChar == "fwn")  // frozen white noise
  {
    long stimulusDuration = atof(serialVals[1]);
    int updateTime = atof(serialVals[2]);
    int nReps = atof(serialVals[3]);
    int randSeedNum = atof(serialVals[4]);
    //Serial.println("Stim: White noise");
    //Serial.flush();
    //  Serial.print("Stim duration: ");
    //  Serial.println(stimulusDuration);
    //  Serial.flush();
    //  Serial.print("Update time: ");
    //   Serial.println(updateTime);
    //  Serial.flush();

    frozenWhiteNoise(updateTime, stimulusDuration, nReps, randSeedNum);

  } else if (FirstChar == "cs")  // contrast switching white noise
  {
    int updateTime = atof(serialVals[1]);
    int switchTime = atof(serialVals[2]);
    int nReps = atof(serialVals[3]);
    float meanVal1 = atof(serialVals[4]);
    float contrastVal1 = atof(serialVals[5]);
    float meanVal2 = atof(serialVals[6]);
    float contrastVal2 = atof(serialVals[7]);
    //Serial.println("Stim: White noise");
    //Serial.flush();
    //  Serial.print("Stim duration: ");
    //  Serial.println(stimulusDuration);
    //  Serial.flush();
    //  Serial.print("Update time: ");
    //   Serial.println(updateTime);
    //  Serial.flush();

    SwitchingWhiteNoise(updateTime, switchTime, nReps, meanVal1, contrastVal1, meanVal2, contrastVal2);
     }
    else if (FirstChar == "se")  // sinusoidal flicker with contrast envelope
    {
     long stimulusDuration = atof(serialVals[1]);
     float frequency = atof(serialVals[2]);
     float envFrequency = atof(serialVals[3]);
    Serial.println("Stim: Sinusoidal env");
    Serial.flush();
    Serial.print("Stim duration: ");
    Serial.println(stimulusDuration);
    Serial.flush();
    Serial.print("Frequency: ");
    Serial.println(frequency);
    Serial.flush();
    Serial.print("Envelope freq: ");
    Serial.println(envFrequency);
    Serial.flush();

    SineContrastConv(stimulusDuration, frequency, envFrequency);

  } else if (FirstChar == "fs")  // frequencySweep stimulus
  {
    float fmin = atof(serialVals[1]);
    float fmax = atof(serialVals[2]);
    float sweepFactorPerSec = atof(serialVals[3]);
    float phaseA = atof(serialVals[4]);
    float phaseB = atof(serialVals[5]);
    float contrastA = atof(serialVals[6]);
    float contrastB = atof(serialVals[7]);

    FrequencySweep(fmin, fmax, sweepFactorPerSec,
                   phaseA, phaseB, contrastA, contrastB);
  } else if (FirstChar == "sfs")  // stepped frequency sweep
  {
    float startFreq = atof(serialVals[1]);
    float endFreq = atof(serialVals[2]);
    float stepFreq = atof(serialVals[3]);
    int cyclesPerFreq = atoi(serialVals[4]);
    float phaseA = atof(serialVals[5]);
    float phaseB = atof(serialVals[6]);
    float contrastA = atof(serialVals[7]);
    float contrastB = atof(serialVals[8]);

    SteppedFrequencySweep(startFreq, endFreq, stepFreq, cyclesPerFreq, phaseA, phaseB, contrastA, contrastB);


  } else if (FirstChar == "sd")  // Set duty cycle of both channels
  {
    float dutyCycle_A = atof(serialVals[1]);
    float dutyCycle_B = atof(serialVals[2]);

    setDutyCycle(dutyCycle_A, dutyCycle_B, TopLumi);

  } else if (FirstChar == "sdt")  // Set duty cycle of both channels for a time period
  {
    float dutyCycle_A = atof(serialVals[1]);
    float dutyCycle_B = atof(serialVals[2]);
    long stimulusDuration = atof(serialVals[3]);

    setDutyCycleTime(dutyCycle_A, dutyCycle_B, stimulusDuration, TopLumi);

  } else if (FirstChar == "gc")  // do gamma correction routine
  {
    float stepSize = atof(serialVals[1]);
    long waitTime = atof(serialVals[2]);
    int nReps = atof(serialVals[3]);
    cycleDutyCycles(stepSize, waitTime, nReps, TopLumi);

  } else if (FirstChar == "useChB")  // apply gamma correction
  {
    useChB = atoi(serialVals[1]);
    if (!useChB) {
      //OCR1B=0;
      Serial.print(F("ChB OFF"));
      Serial.print("\n");
    } else {
      Serial.print(F("ChB ON"));
      Serial.print("\n");
    };

  } else if (FirstChar == "useChA")  // apply gamma correction
  {
    useChA = atoi(serialVals[1]);
    if (!useChA) {
      //OCR1A=0;
      Serial.print(F("ChA OFF"));
      Serial.print("\n");
    } else {
      Serial.print(F("ChA ON"));
      Serial.print("\n");
    };
  }  //else if (FirstChar == "stat") {
  //getStatus();

  //}
  else if (FirstChar == "ana") {
    readAnalogVals();

  } else if (FirstChar == "agc")  // apply gamma correction
  {
    uint8_t lutIndex = atoi(serialVals[1]);
    if (lutIndex == 1) {
      currentChALUT = ChA1LUT;
      currentChBLUT = ChA1LUT;
      Serial.print(F("LUT 1 SELECTED"));
      Serial.print("\n");

    } else if (lutIndex == 2) {
      currentChALUT = ChA2LUT;
      currentChBLUT = ChA2LUT;
      Serial.print(F("LUT 2 SELECTED"));
      Serial.print("\n");
    }
  } else  // not valid stimulus code
  {
    Serial.print(FirstChar);
    Serial.print(F(" is an invalid stimulus code - make sure you are using carriage return line ending"));
    Serial.print("\n");
  }
  //memset('\0', receivedChars, sizeof(receivedChars));
  memset(receivedChars, '\0', sizeof(receivedChars));
  //Serial.print("rc: ");
  //Serial.println(receivedChars);
}


///////////////////////////////////// SINEWAVE FLICKER  //////////////////////////////////////
// Function to generate a sine wave table
void generateSineWaveTable(long TOP) {
  for (int i = 0; i < TABLE_SIZE; i++) {
    // Calculate the sine wave value (scaled between 0 and TOP)
    float angle = (2.0 * PI * i) / TABLE_SIZE;                        // Angle in radians
    sineWaveTable[i] = (uint16_t)((sin(angle) + 1.0) * (TOP / 2.0));  // Scale to 0-TOP, un-corrected sinewave
  }
}

void outputSinewave(float sinewaveFrequency, long duration, float phaseA, float phaseB, float contrastA, float contrastB) {

  tableIndexA = phaseA * (float)TABLE_SIZE;  // Start at the beginning of the sine wave table
  tableIndexB = phaseB * (float)TABLE_SIZE;

  contrastMultA = contrastA;
  contrastMultB = contrastB;

  // first do some calculations to find the update interval and step size for Timer3 interrupts
  // Calculate the PWM cycle time in microseconds
  float pwmCycleTime = (2.0 * TOP) / (float)(CLOCK_FREQ / prescaler);  // Time per PWM cycle in seconds
  pwmCycleTime *= 1e6;                                                 // Convert seconds to microseconds

  // Calculate the base update interval for the sinewave frequency
  float baseUpdateInterval = 1.0 / (sinewaveFrequency * TABLE_SIZE);  // Time per table update in seconds
  baseUpdateInterval *= 1e6;                                          // Convert to microseconds

  // Find the smallest step size that is a factor of TABLE_SIZE
  stepSize = TABLE_SIZE;  // Start with the maximum possible step size
  for (int i = 1; i <= TABLE_SIZE; i++) {
    if ((TABLE_SIZE % i == 0) && (baseUpdateInterval * i >= pwmCycleTime)) {
      stepSize = i;
      break;  // Stop at the first valid (smallest) step size
    }
  }

  // Recalculate the effective update interval based on the step size
  float updateInterval = baseUpdateInterval * stepSize;
  float updateFrequency = 1e6 / updateInterval;  // update frequency for timer3 interrupt

  configureTimer3Interrupt(updateFrequency);  // Configure timer3 interrupt to updateFrequency
  PORTC |= (1 << PORTC6);                     // set Pin 5 HIGH
  PORTD |= (1 << PIND4);                      // set Pin 4 HIGH

  long startTime = millis();  // Record the start time
  // set timer3 interrupt callback function to play the sinewave
  setTimer3Callback(sinewaveInterrupt);

  // Loop until the specified duration has elapsed
  while (millis() - startTime < duration) {
    delayMicroseconds(1);  //wait for time to end
  }

  stopTimer3Interrupt();    // finish playing sinewave
  PORTD &= ~(1 << PIND4);   // Ensure Pin 4 is set to LOW by changing register directly
  PORTC &= ~(1 << PORTC6);  // Ensure Pin 5 is set to LOW
  Serial.print("-1");
  Serial.print("\n");

  if (useChA) { setChA(TopLumi / 2); }  // Set pin 9 to 50% duty cycle as default
  if (useChB) { setChB(TopLumi / 2); }  // Set pin 10 to 50% duty cycle as default
}

// sinewave interrupt function
void sinewaveInterrupt() {
  //static int tableIndexA = 0;  // Start at the beginning of the sine wave table
  // Update PWM duty cycle with the next sine wave value

  //ocrVal = MidLumi + ((sineWaveTable[tableIndexA] - MidLumi) * (contrastA));
  float ocrValA = MidLumi + ((sineWaveTable[tableIndexA] - MidLumi) * (contrastMultA));
  float ocrValB = MidLumi + ((sineWaveTable[tableIndexB] - MidLumi) * (contrastMultB));


  if (useChA) { setChA(ocrValA); }  //
  if (useChB) { setChB(ocrValB); }  //
  //Serial.print(OCR1A);
  //Serial.print(',');

  // Update the table index (wrap around if necessary)
  tableIndexA = tableIndexA + stepSize;
  if (tableIndexA >= TABLE_SIZE) {
    PORTD ^= (1 << PIND4);      // Toggle Pin 4 if sine wave cycle finished
    tableIndexA -= TABLE_SIZE;  // wrap table
    completedCycles++;          // increment completed cycles counter
  }

  // Update the table index (wrap around if necessary)
  tableIndexB = tableIndexB + stepSize;
  if (tableIndexB >= TABLE_SIZE) {
    tableIndexB -= TABLE_SIZE;  // wrap table
  }
}

/////////////////////////////////// SINE WAVE FLICKER WITH CONTRAST ENVELOPE //////////////////////////
void SineContrastConv(float duration, float sinewaveFrequency, float envelopeFreq) {

// first do some calculations to find the update interval and step size for Timer3 interrupts
// Calculate the PWM cycle time in microseconds
  float pwmCycleTime = (2.0 * TOP) / (float)(CLOCK_FREQ / prescaler);  // Time per PWM cycle in seconds
  pwmCycleTime *= 1e6;                                                 // Convert seconds to microseconds

// Calculate the base update interval for the sinewave frequency
 float baseUpdateInterval = 1.0 / (sinewaveFrequency * TABLE_SIZE);  // Time per table update in seconds
 baseUpdateInterval *= 1e6;                                          // Convert to microseconds

// Find the smallest step size that is a factor of TABLE_SIZE
 stepSize = TABLE_SIZE;  // Start with the maximum possible step size
  for (int i = 1; i <= TABLE_SIZE; i++) {
    if ((TABLE_SIZE % i == 0) && (baseUpdateInterval * i >= pwmCycleTime)) {
      stepSize = i;
      break;  // Stop at the first valid (smallest) step size
    }
  }

  tableIndexA = 0;      // start at beginning of sinewave table
 tableEnvIndex = 191;  // start at 0 contrast
 envCount = 0;         // contrast envelope counter
 contrastMult = 0;
 nEnvCounts = 0;


// Recalculate the effective update interval based on the step size
 float updateInterval = baseUpdateInterval * stepSize;
//Serial.print("req update: ");
//Serial.println(updateInterval);
 float updateFrequency = 1e6 / updateInterval;  // update frequency for timer3 interrupt
 nEnvCounts = (int)(sinewaveFrequency / envelopeFreq);

// Configure timer3 interrupt to updateFrequency
 configureTimer3Interrupt(updateFrequency);
 PORTC |= (1 << PORTC6);     // Stim on pin 5
 long startTime = millis();  // Record the start time
// set timer3 interrupt callback function to play the sinewave
  setTimer3Callback(sinewaveEnvelopeInterrupt);

// Loop until the specified duration has elapsed
 while (millis() - startTime < duration) {
    delayMicroseconds(1);  //wait for time to end
  }
  stopTimer3Interrupt();    // finish playing sinewave
 PORTD &= ~(1 << PIND4);   // Ensure Pin 4 is set to LOW by changing register directly
 PORTC &= ~(1 << PORTC6);  // Ensure Pin 5 is set to LOW
 Serial.println("-1");
 Serial.flush();
 if (useChA) { setChA(TopLumi / 2); }  // Set pin 9 to 50% duty cycle as default
 if (useChB) { setChB(TopLumi / 2); }  // Set pin 10 to 50% duty cycle as default}
}


// sinewave contrast envelope interrupt function
void sinewaveEnvelopeInterrupt() {
unsigned long startTime = micros();
// Update PWM duty cycle with the next sine wave value
  uint16_t ocrVal = MidLumi + ((sineWaveTable[tableIndexA] - MidLumi) * (contrastMult));

 if (useChA) { setChA(ocrVal); }  //
  if (useChB) { setChB(ocrVal); }  //

// update sinewave table index based on interrupt frequency
 tableIndexA = tableIndexA + stepSize;
 if (tableIndexA >= TABLE_SIZE) tableIndexA -= TABLE_SIZE;  // wrap table index


 //update counter for contrast envelope
 envCount = envCount + 1;
 if (envCount > nEnvCounts - 1)  // check if time to update contrast value
 {
  envCount = 0;                       // reset
   tableEnvIndex = tableEnvIndex + 1;  // incremenet contrast LUT index

   if (tableEnvIndex >= TABLE_SIZE) {
    PORTD ^= (1 << PIND4);        // Toggle Pin 4 if envelope cycle finished
     tableEnvIndex -= TABLE_SIZE;  // wrap tableEnvIndex
  }
   contrastMult = sineWaveTable[tableEnvIndex] / float(TopLumi);  // update contrast multiplier
 }
unsigned long duration = micros() - startTime; // Measure execution time
//Serial.println(duration); // Print execution time
}

///////////////////////////////////  FREQUENCY SWEEP FUNCTIONS//////////////////////////////////

void SteppedFrequencySweep(float startFreq, float endFreq, float stepFreq, int cyclesPerFreq, float phaseA, float phaseB, float contrastA, float contrastB) {

  // Set initial phase and contrast
  tableIndexA = (int)(phaseA * (float)TABLE_SIZE) % TABLE_SIZE;
  tableIndexB = (int)(phaseB * (float)TABLE_SIZE) % TABLE_SIZE;
  contrastMultA = contrastA;
  contrastMultB = contrastB;

  PORTC |= (1 << PORTC6);  // set Pin 5 HIGH
  PORTD |= (1 << PIND4);   // set Pin 4 HIGH

  float currentFreq = startFreq;
  bool sweepingUp = startFreq <= endFreq;
  stepFreq = abs(stepFreq); 

  // Set the callback function once
  setTimer3Callback(sinewaveInterrupt);

  while ((sweepingUp && currentFreq <= endFreq) || (!sweepingUp && currentFreq >= endFreq)) {
    
    // 1. Calculate the timer math for the current frequency
    float pwmCycleTime = (2.0 * TOP) / (float)(CLOCK_FREQ / prescaler) * 1e6;
    float baseUpdateInterval = (1.0 / (currentFreq * TABLE_SIZE)) * 1e6;

    stepSize = TABLE_SIZE;
    for (int i = 1; i <= TABLE_SIZE; i++) {
      if ((TABLE_SIZE % i == 0) && (baseUpdateInterval * i >= pwmCycleTime)) {
        stepSize = i;
        break;
      }
    }

    float updateInterval = baseUpdateInterval * stepSize;
    float updateFrequency = 1e6 / updateInterval;

    // 2. Reset cycle counter BEFORE configuring the new timer speed
    completedCycles = 0; 
    
    // 3. Apply the new frequency (changes how fast the interrupt fires)
    configureTimer3Interrupt(updateFrequency);

    // 4. Wait for EXACTLY 'cyclesPerFreq' wrap-arounds
    while (completedCycles < cyclesPerFreq) {
      // The microcontroller just hangs out here.
      // All the precision work is happening in the background via Timer 3.
      delayMicroseconds(1); 
    }

    // 5. Step the frequency for the next loop
    if (sweepingUp) {
      currentFreq += stepFreq;
    } else {
      currentFreq -= stepFreq;
    }
  }

  // Cleanup after sweep is done
  stopTimer3Interrupt();
  PORTD &= ~(1 << PIND4);   // Ensure Pin 4 is LOW
  PORTC &= ~(1 << PORTC6);  // Ensure Pin 5 is LOW
  Serial.print("-1\n");

  if (useChA) { setChA(TopLumi / 2); } 
  if (useChB) { setChB(TopLumi / 2); } 
}

/**
 * Outputs a sine wave with an approximated exponential frequency sweep for one full cycle.
 * The frequency sweeps from fmin to fmax, then back to fmin, using simplified step multiplication.
 * The total duration is calculated based on fmin, fmax, and sweepFactorPerSec using log() once.
 *
 */
void FrequencySweep(float fmin, float fmax, float sweepFactorPerSec,
                    float phaseA, float phaseB, float contrastA, float contrastB) {

  // Ensure fmin and fmax are positive
  if (fmin <= 0.0f) fmin = 0.001f;
  if (fmax <= 0.0f) fmax = 0.001f;

  if (fmin > fmax) {
    float temp = fmin;
    fmin = fmax;
    fmax = temp;
  }

  tableIndexA = (int)(phaseA * (float)TABLE_SIZE) % TABLE_SIZE;
  tableIndexB = (int)(phaseB * (float)TABLE_SIZE) % TABLE_SIZE;
  contrastMultA = contrastA;
  contrastMultB = contrastB;

  float pwmCycleTimeUs = (2.0f * TOP) / (float)(CLOCK_FREQ / prescaler);
  pwmCycleTimeUs *= 1e6f;

  unsigned long totalDurationUs = 0;
  unsigned long timeForOneWaySweepUs = 0;
  bool isSweeping = false;

  // Calculate dt in seconds for frequency update steps
  const float dt_sec_step = (float)TARGET_RECONFIG_INTERVAL_US / 1000000.0f;
  float step_mult_change_factor = 0.0f;  // This is M * dt

  // Check for valid sweep conditions. sweepFactorPerSec is M from f(t) = f0 * exp(M*t)
  // Time to sweep from fmin to fmax is t = log(fmax/fmin) / M
  if (sweepFactorPerSec > 0.000001f && fmax > fmin && (fmax / fmin) > 1.000001f) {
    float timeForOneWaySweepSec_calc = log(fmax / fmin) / sweepFactorPerSec;

    if (timeForOneWaySweepSec_calc > 0.0000001f) {
      isSweeping = true;
      timeForOneWaySweepUs = (unsigned long)(timeForOneWaySweepSec_calc * 1000000.0f);
      totalDurationUs = 2 * timeForOneWaySweepUs;
      step_mult_change_factor = sweepFactorPerSec * dt_sec_step;  // M*dt

      if (totalDurationUs == 0 && timeForOneWaySweepSec_calc > 0.0000001f) {
        totalDurationUs = 2;  // Ensure at least minimal duration if calculated sweep is very short but valid
        if (timeForOneWaySweepUs == 0) timeForOneWaySweepUs = 1;
      } else if (totalDurationUs == 0) {
        isSweeping = false;
      }
    } else {
      isSweeping = false;  // Calculated time is too short or invalid
    }
  }

  if (!isSweeping) {
    totalDurationUs = TARGET_RECONFIG_INTERVAL_US;
    timeForOneWaySweepUs = totalDurationUs / 2;  // Not strictly used but avoids being zero
                                                 // For fixed frequency, step_mult_change_factor remains 0, so freq won't change.
  }

  if (totalDurationUs == 0) {
    totalDurationUs = TARGET_RECONFIG_INTERVAL_US;
    timeForOneWaySweepUs = totalDurationUs / 2;
    isSweeping = false;
  }

  float actual_calc_freq = fmin;       // This holds the frequency to be set for the current/next interval
  bool sweeping_up = true;             // State variable to track sweep direction (more robustly handled by time check)
  bool just_switched_to_down = false;  // Flag to handle the exact moment of switching to downward sweep

  // get initial timer interrupt frequency
  float currentSinewaveFrequency = actual_calc_freq;

  // Clamp and ensure positive frequency before using it for calculations
  if (currentSinewaveFrequency > fmax) currentSinewaveFrequency = fmax;
  if (currentSinewaveFrequency < fmin) currentSinewaveFrequency = fmin;
  if (currentSinewaveFrequency <= 0.0f) currentSinewaveFrequency = 0.0001f;


  // --- Recalculate timer parameters based on currentSinewaveFrequency ---
  float sinewaveFreqToUseForCalc = currentSinewaveFrequency;
  float baseUpdateIntervalUs = (1.0f / (sinewaveFreqToUseForCalc * TABLE_SIZE)) * 1e6f;
  int calculatedStepSize = TABLE_SIZE;
  for (int i = 1; i <= TABLE_SIZE; i++) {
    if ((TABLE_SIZE % i == 0) && (baseUpdateIntervalUs * i >= pwmCycleTimeUs)) {
      calculatedStepSize = i;
      break;
    }
  }

  //noInterrupts();
  stepSize = calculatedStepSize;
  //interrupts();

  float effectiveTimerUpdateIntervalUs = baseUpdateIntervalUs * calculatedStepSize;
  float timerInterruptFrequencyHz;
  timerInterruptFrequencyHz = 1e6f / effectiveTimerUpdateIntervalUs;
  configureTimer3Interrupt(timerInterruptFrequencyHz);

  PORTC |= (1 << PORTC6);  // set Pin 5 HIGH
  PORTD |= (1 << PIND4);   // set pin 4 high initially, toggles every sinewave cycle
  unsigned long startTimeUs = micros();
  setTimer3Callback(sinewaveInterrupt);

  // Main loop
  while (micros() - startTimeUs < totalDurationUs) {
    unsigned long loopIterationStartTimeUs = micros();

    // currentSinewaveFrequency is the frequency for the current Timer3 configuration
    float currentSinewaveFrequency = actual_calc_freq;

    // Clamp and ensure positive frequency before using it for calculations
    if (currentSinewaveFrequency > fmax) currentSinewaveFrequency = fmax;
    if (currentSinewaveFrequency < fmin) currentSinewaveFrequency = fmin;
    if (currentSinewaveFrequency <= 0.0f) currentSinewaveFrequency = 0.0001f;


    // --- Recalculate timer parameters based on currentSinewaveFrequency ---
    float sinewaveFreqToUseForCalc = currentSinewaveFrequency;
    float baseUpdateIntervalUs = (1.0f / (sinewaveFreqToUseForCalc * TABLE_SIZE)) * 1e6f;
    int calculatedStepSize = TABLE_SIZE;
    for (int i = 1; i <= TABLE_SIZE; i++) {
      if ((TABLE_SIZE % i == 0) && (baseUpdateIntervalUs * i >= pwmCycleTimeUs)) {
        calculatedStepSize = i;
        break;
      }
    }

    //noInterrupts();
    stepSize = calculatedStepSize;
    //interrupts();

    float effectiveTimerUpdateIntervalUs = baseUpdateIntervalUs * calculatedStepSize;
    float timerInterruptFrequencyHz;
    if (effectiveTimerUpdateIntervalUs <= 0.000001f) {
      timerInterruptFrequencyHz = 0;
    } else {
      timerInterruptFrequencyHz = 1e6f / effectiveTimerUpdateIntervalUs;
    }

    if (timerInterruptFrequencyHz > 0) {
      configureTimer3Interrupt(timerInterruptFrequencyHz);
    } else {
      stopTimer3Interrupt();
    }

    // --- Update actual_calc_freq for the NEXT interval ---
    if (isSweeping) {
      unsigned long elapsedTimeTotalUs = loopIterationStartTimeUs - startTimeUs;

      if (elapsedTimeTotalUs < timeForOneWaySweepUs) {  // Sweeping up phase
        actual_calc_freq = actual_calc_freq * (1.0f + step_mult_change_factor);
        if (actual_calc_freq > fmax) {
          actual_calc_freq = fmax;
        }
        just_switched_to_down = true;  // Reset flag, ready for when down sweep starts
      } else {                         // Sweeping down phase
        if (just_switched_to_down) {   // First time entering down sweep based on time
          actual_calc_freq = fmax;     // Ensure we start down sweep precisely from fmax
          just_switched_to_down = false;
        }
        actual_calc_freq = actual_calc_freq * (1.0f - step_mult_change_factor);
        if (actual_calc_freq < fmin) {
          actual_calc_freq = fmin;
        }
      }
      // Final safety clamp for actual_calc_freq after multiplication
      if (actual_calc_freq <= 0.0f) actual_calc_freq = fmin;                                                    // Prevent zero/negative
      else if (actual_calc_freq > fmax && elapsedTimeTotalUs >= timeForOneWaySweepUs) actual_calc_freq = fmax;  // Cap if overshot fmax during down phase start
      else if (actual_calc_freq < fmin) actual_calc_freq = fmin;


    } else {
      actual_calc_freq = fmin;  // Fixed frequency, no change
    }

    // --- Loop Pacing ---
    unsigned long loopProcessingTimeUs = micros() - loopIterationStartTimeUs;
    if (loopProcessingTimeUs < TARGET_RECONFIG_INTERVAL_US) {
      long delayNeededUs = TARGET_RECONFIG_INTERVAL_US - loopProcessingTimeUs;
      unsigned long currentTimeInLoopUs = micros() - startTimeUs;
      if (currentTimeInLoopUs + delayNeededUs > totalDurationUs) {
        if (totalDurationUs > currentTimeInLoopUs) {
          delayNeededUs = totalDurationUs - currentTimeInLoopUs;
        } else {
          delayNeededUs = 0;
        }
      }
      if (delayNeededUs > 0) {
        delayMicroseconds(delayNeededUs);
      }
    }
  }

  stopTimer3Interrupt();    // finish playing sinewave
  PORTD &= ~(1 << PIND4);   // Ensure Pin 4 is set to LOW by changing register directly
  PORTC &= ~(1 << PORTC6);  // Ensure Pin 5 is set to LOW
  Serial.print("-1");
  Serial.print("\n");

  if (useChA) { setChA(TopLumi / 2); }  // Set pin 9 to 50% duty cycle as default
  if (useChB) { setChB(TopLumi / 2); }  // Set pin 10 to 50% duty cycle as default
}




/////////////////////////////////// WHITE NOISE PWM FUNCTIONS //////////////////////////////

// --- Function to generate N(0,1) using CLT ---
float generateGaussianCLT() {
  const int N = CLT_N;
  long sum_random = 0;
  for (int i = 0; i < N; i++) {
    // sum_random += random(RANDOM_UPPER_BOUND); // Original
    sum_random += fast_rand32() & (RANDOM_UPPER_BOUND - 1);  // Faster core PRNG (RANDOM_UPPER_BOUND must be power of 2, e.g. 1024)
  }
  float gaussian_approx_zero_mean_unit_variance =
    ((float)sum_random - SUM_EXPECTED_MEAN) * NORMALIZE_SCALE_FACTOR;
  return gaussian_approx_zero_mean_unit_variance;
}

void whiteNoise(long updateTime, long duration, float frac_target_mean, float frac_target_std) {

  //printSequenceNum=1; // for serial debugging
  float updateFrequency = 1e3 / updateTime;
  configureTimer3Interrupt(updateFrequency);

  target_mean = TopLumi * frac_target_mean;
  target_std = TopLumi * frac_target_std;

  // clamping limits
  map_float_min = (float)0.0;
  map_float_max = (float)TopLumi;

  Serial.print("TOP: ");
  Serial.print(TopLumi);
  Serial.print("\n");

  PORTC |= (1 << PORTC6);     // Stim on pin 5
  long startTime = millis();  // Record the start time

  // set timer3 interrupt callback function to play the sinewave
  setTimer3Callback(whiteNoiseInterrupt);

  // Loop until the specified duration has elapsed
  while (millis() - startTime < duration) {
    delayMicroseconds(1);  //wait for time to end
  }

  stopTimer3Interrupt();  // stop white noise
  delay(updateTime);
  PORTD &= ~(1 << PIND4);   // Ensure Pin 4 is set to LOW by changing register directly
  PORTC &= ~(1 << PORTC6);  // Ensure Pin 5 is set to LOW
  Serial.print("-1");
  Serial.print("\n");
  Serial.flush();
  if (useChA) { setChA(TopLumi / 2); }  // Set pin 9 to 50% duty cycle as default
  if (useChB) { setChB(TopLumi / 2); }  // Set pin 10 to 50% duty cycle as default
}

void SwitchingWhiteNoise(long updateTime, unsigned long switchTime, int nReps, float meanVal1, float contrastVal1, float meanVal2, float contrastVal2) {

  float updateFrequency = 1e3 / updateTime;
  configureTimer3Interrupt(updateFrequency);

  long duration = switchTime * 2 * nReps;
  Serial.print("switch time:");
  Serial.print(switchTime);
  Serial.print("\n");
  Serial.print("duration:");
  Serial.print(duration);
  Serial.print("\n");

  // start with dist1 values
  target_mean = TopLumi * meanVal1;
  target_std = TopLumi * contrastVal1;
  currentDist = 1;

  // clamping limits
  map_float_min = (float)0.0;
  map_float_max = (float)TopLumi;

  Serial.print("TOP: ");
  Serial.print(TopLumi);
  Serial.print("\n");

  PORTC |= (1 << PORTC6);  // Stim on pin 5

  long startTime = millis();  // Record the start time

  // set timer3 interrupt callback function to play the sinewave
  setTimer3Callback(whiteNoiseInterrupt);

  unsigned long switchPreviousMillis = millis();  // will store last time distribution switched
  unsigned long switchCurrentMillis = millis();

  // Loop until the specified duration has elapsed
  while (millis() - startTime < duration) {
    delayMicroseconds(1);                                            //wait for time to end
    switchCurrentMillis = millis();                                  // get current time to check for whether to change distirbution
    if (switchCurrentMillis - switchPreviousMillis >= switchTime) {  // if time to switch distributions
                                                                     //Serial.println(switchCurrentMillis);

      switchPreviousMillis = switchCurrentMillis;  // reset switchPreviousMillis

      if (currentDist == 1) {
        target_mean = TopLumi * meanVal2;
        target_std = TopLumi * contrastVal2;
        currentDist = 2;
        //Serial.println(target_std);
      } else if (currentDist == 2) {
        target_mean = TopLumi * meanVal1;
        target_std = TopLumi * contrastVal1;
        currentDist = 1;
        //Serial.println(target_std);
      }

      PORTC ^= (1 << PORTC6);  //Toggle Pin 5
    }
  }

  stopTimer3Interrupt();  // stop white noise
  delay(updateTime);
  PORTD &= ~(1 << PIND4);   // Ensure Pin 4 is set to LOW by changing register directly
  PORTC &= ~(1 << PORTC6);  // Ensure Pin 5 is set to LOW
  Serial.print("-1");
  Serial.print("\n");
  Serial.flush();
  if (useChA) { setChA(TopLumi / 2); }  // Set pin 9 to 50% duty cycle as default
  if (useChB) { setChB(TopLumi / 2); }  // Set pin 10 to 50% duty cycle as default
}

// whitenoise interrupt function
void whiteNoiseInterrupt() {

  // Update PWM duty cycle with the next random value
  if (useChA) { setChA(finalRandNumber_A); }
  if (useChB) { setChB(finalRandNumber_B); }

  if (useChA) {
    float zA = generateGaussianCLT();
    float float_gaussian_val_A = target_mean + target_std * zA;
    if (float_gaussian_val_A < map_float_min) { float_gaussian_val_A = map_float_min; }
    if (float_gaussian_val_A > map_float_max) { float_gaussian_val_A = map_float_max; }
    finalRandNumber_A = (uint16_t)float_gaussian_val_A;
  }

  if (useChB) {
    float zB = generateGaussianCLT();
    float float_gaussian_val_B = target_mean + target_std * zB;
    if (float_gaussian_val_B < map_float_min) { float_gaussian_val_B = map_float_min; }
    if (float_gaussian_val_B > map_float_max) { float_gaussian_val_B = map_float_max; }
    finalRandNumber_B = (uint16_t)float_gaussian_val_B;
  }

  PIND = (1 << PIND4);  // alternate PIN 4 value indicator pin
  Serial.print(printSequenceNum);
  Serial.print(",");
  Serial.print(finalRandNumber_A);
  Serial.print(",");
  Serial.print(finalRandNumber_B);
  Serial.print("\n");
  //Serial.print(",");

  printSequenceNum++;
}


/////////////////////////////////// FROZEN WHITE NOISE PWM FUNCTIONS //////////////////////////////

void frozenWhiteNoise(int updateTime, long duration, long nReps, int randSeedNum) {

  tableIndexFWN = 0;  // start at beginning of frozen white noise segment
  float updateFrequency = 1e3 / updateTime;
  configureTimer3Interrupt(updateFrequency);

  long totalDuration = duration * nReps;
  Serial.print("LD: ");
  Serial.print(totalDuration);
  Serial.print("\n");
  xorshift32_state = randSeedNum;  // for reproducible random sequence across different stimulus blocks


  float frac_target_mean = 0.5;
  float frac_target_std = 0.2;
  target_mean = TopLumi * frac_target_mean;
  target_std = TopLumi * frac_target_std;

  // clamping limits
  map_float_min = (float)0.0;
  map_float_max = (float)TopLumi;

  for (int i = 0; i < FWN_TABLE_SIZE; i++) {

    float zA = generateGaussianCLT();
    float float_gaussian_val_A = target_mean + target_std * zA;
    if (float_gaussian_val_A < map_float_min) { float_gaussian_val_A = map_float_min; }
    if (float_gaussian_val_A > map_float_max) { float_gaussian_val_A = map_float_max; }
    finalRandNumber_A = (uint16_t)float_gaussian_val_A;
    frozenWhiteNoiseTable[i] = finalRandNumber_A;
  }

  Serial.print("TOP: ");
  Serial.print(TopLumi);
  Serial.print("\n");
  Serial.flush();

  long startTime = millis();  // Record the start time
  // set timer3 interrupt callback function to play the sinewave
  PORTC |= (1 << PORTC6);  // Stim on pin 5
  setTimer3Callback(frozenWhiteNoiseInterrupt);

  // Loop until the specified duration has elapsed
  while (millis() - startTime < totalDuration) {
    delayMicroseconds(1);  //wait for time to end
  }

  stopTimer3Interrupt();  // stop white noise
  delay(updateTime);
  PORTD &= ~(1 << PIND4);   // Ensure Pin 4 is set to LOW by changing register directly
  PORTC &= ~(1 << PORTC6);  // Ensure Pin 5 is set to LOW
  Serial.print("-1");
  Serial.print("\n");
  Serial.flush();
  if (useChA) { setChA(TopLumi / 2); }  // Set pin 9 to 50% duty cycle as default
  if (useChB) { setChB(TopLumi / 2); }  // Set pin 10 to 50% duty cycle as default
}

// whitenoise interrupt function
void frozenWhiteNoiseInterrupt() {

  // Update PWM duty cycle with the frozen white noise value
  if (useChA) { setChA(frozenWhiteNoiseTable[tableIndexFWN]); }
  if (useChB) { setChB(frozenWhiteNoiseTable[tableIndexFWN]); }
  PIND = (1 << PIND4);  // alternate PIN 4 value indicator pin


  //Serial.print("ti: ");
  //Serial.println(tableIndexA);
  Serial.print(frozenWhiteNoiseTable[tableIndexFWN]);
  //Serial.print(",");
  Serial.print("\n");
  Serial.flush();

  // Update the table index (wrap around at actual white noise table size)
  tableIndexFWN = tableIndexFWN + 1;  // increment table index
  if (tableIndexFWN >= FWN_TABLE_SIZE) {
    PORTC ^= (1 << PORTC6);           //Toggle Pin 5 when table finishes
    tableIndexFWN -= FWN_TABLE_SIZE;  // wrap table
  }
}

/////////////////////////////////// SOME GENERIC PWM FUNCTIONS ///////////////////////////////////////////

// set gamma-corrected output of channel A (pin 9)
void setChA(uint16_t ocrValue) {
  OCR1A = pgm_read_word_near(currentChALUT + ocrValue);
}

// set gamma-corrected output of channel B (pin 10)
void setChB(uint16_t ocrValue) {
  OCR1B = pgm_read_word_near(currentChBLUT + ocrValue);
}

// function to artifically lower the max PWM duty cycle. (i.e. TopMultiplier=0.5 means max duty cycle of 50%)
// other functions will work as normal but scale to this TOP value
//void SetTopLumi(float TopMultiplier) {

//if (TopMultiplier > 1) {
//    TopMultiplier = 1;
//  } else if (TopMultiplier <= 0) {
//    TopMultiplier = 1;
//  }

// get new TOP value to use
//  TopLumi = float(TOP) * TopMultiplier;
//  MidLumi = TopLumi / 2;

// Generate the sine wave LUT based on the Timer1 config
// generateSineWaveTable(TopLumi);

// initialise random number to 50% duty cyle
// randNumber = TopLumi / 2;

// Set pin 9 to 50% duty cycle as default
// if (useChA) { setChA(TopLumi / 2); }  // Set pin 9 to 50% duty cycle as default
// if (useChB) { setChB(TopLumi / 2); }  // Set pin 10 to 50% duty cycle as default
//}


// set the duty cycle manually until a new value is requested
void setDutyCycle(float dutyCyclePercentage_A, float dutyCyclePercentage_B, long TopLumi) {
  // Constrain the duty cycle percentage between 0% and 100%
  if (dutyCyclePercentage_A < 0.0) dutyCyclePercentage_A = 0.0;
  if (dutyCyclePercentage_A > 100.0) dutyCyclePercentage_A = 100.0;
  // Calculate the OCR1 value based on the duty cycle and TOP
  uint16_t ocrValueA = (long)((dutyCyclePercentage_A / 100.0) * TopLumi);

  if (dutyCyclePercentage_B < 0.0) { dutyCyclePercentage_B = 0.0; };
  if (dutyCyclePercentage_B > 100.0) { dutyCyclePercentage_B = 100.0; }
  // Calculate the OCR1 value based on the duty cycle and TOP
  uint16_t ocrValueB = (long)((dutyCyclePercentage_B / 100.0) * TopLumi);

  // Set OCR1A to control the duty cycle
  PORTD ^= (1 << PIND4);  // toggle pin 4 whenever duty cycles are changed
  //if (useChA)
  setChA(ocrValueA);
  //if (useChB)
  setChB(ocrValueB);
  //Serial.print(ocrValue);
  //Serial.print(',');
  //Serial.println(OCR1A);
}


// set the duty cycle manually until a new value is requested
void setDutyCycleTime(float dutyCyclePercentage_A, float dutyCyclePercentage_B, long duration, long TopLumi) {
  // Constrain the duty cycle percentage between 0% and 100%
  if (dutyCyclePercentage_A < 0.0) dutyCyclePercentage_A = 0.0;
  if (dutyCyclePercentage_A > 100.0) dutyCyclePercentage_A = 100.0;
  // Calculate the OCR1 value based on the duty cycle and TOP
  uint16_t ocrValueA = (long)((dutyCyclePercentage_A / 100.0) * TopLumi);

  if (dutyCyclePercentage_B < 0.0) { dutyCyclePercentage_B = 0.0; };
  if (dutyCyclePercentage_B > 100.0) { dutyCyclePercentage_B = 100.0; }
  // Calculate the OCR1 value based on the duty cycle and TOP
  uint16_t ocrValueB = (long)((dutyCyclePercentage_B / 100.0) * TopLumi);


  long startTime = millis();          // Record the start time
  PORTC |= (1 << PORTC6);             // Stim on pin 5
  PORTD |= (1 << PIND4);              // match pin 5 on pin 4 for this stimulus
  if (useChA) { setChA(ocrValueA); }  // Set pin 9 to 50% duty cycle as default
  if (useChB) { setChB(ocrValueB); }  // Set pin 10 to 50% duty cycle as default
  // Loop until the specified duration has elapsed
  while (millis() - startTime < duration) {
    delayMicroseconds(1);  //wait for time to end
  }
  PORTD &= ~(1 << PIND4);   // Ensure Pin 4 is set to LOW by changing register directly
  PORTC &= ~(1 << PORTC6);  // Ensure Pin 5 is set to LOW
  Serial.print("-1");
  Serial.print("\n");
  Serial.flush();
  if (useChA) { setChA(TopLumi / 2); }  // Set pin 9 to 50% duty cycle as default
  if (useChB) { setChB(TopLumi / 2); }  // Set pin 10 to 50% duty cycle as default
}


// Run through duty cycles to perform gamma correction
void cycleDutyCycles(float stepSize, float waitTime, int nReps, long TopLumi) {
  float dutyCycle = 0;

  for (int irep = 0; irep < nReps; irep++) {
    while (dutyCycle <= 1) {
      Serial.print(dutyCycle);
      Serial.print("\n");
      long ocrValue = (long)(dutyCycle * TopLumi);
      //OCR1A = ocrValue;
      if (useChA) { setChA(ocrValue); }  // Set pin 9 to 50% duty cycle as default
      if (useChB) { setChB(ocrValue); }  // Set pin 10 to 50% duty cycle as default    delay(waitTime);
      dutyCycle = dutyCycle + stepSize;
      delay(waitTime);
    }
  }
  Serial.print("-1");
  Serial.print("\n");
  // Set pin 9 to 50% duty cycle as default
  if (useChA) { setChA(TopLumi / 2); }  // Set pin 9 to 50% duty cycle as default
  if (useChB) { setChB(TopLumi / 2); }  // Set pin 10 to 50% duty cycle as default
}

//void getStatus() {
//  Serial.print(F("PWM FREQ: "));
///  Serial.println(desiredPWMFrequency);
/// Serial.print(F("TOP: "));
//  Serial.println(TOP);
// Serial.print(F("TopLumi: "));
// // Serial.println(TopLumi);
//  Serial.print(F("Duty cycle: "));
///  Serial.println(dutyCycle);
//  Serial.print(F("Current OCR1A: "));
//  Serial.println(OCR1A);
//  Serial.print(F("Current OCR1B: "));
//  Serial.println(OCR1B);
//}


void readAnalogVals() {
  //setDutyCycle(100, 100, TopLumi);  //set max duty cycle to get clean readings
  bool keepReading = true;
  const unsigned long interval = 100;  // 100ms interval
  unsigned long previousMillis = millis();

  int analogValue0 = analogRead(A0);
  int analogValue1 = analogRead(A1);

  while (keepReading) {
    // Check if the interval has passed
    unsigned long currentMillis = millis();
    if (currentMillis - previousMillis >= interval) {
      previousMillis = currentMillis;

      // Read and print analog values
      analogValue0 = analogRead(A0);
      analogValue1 = analogRead(A1);
      Serial.print(analogValue0);
      Serial.print(",");
      Serial.println(analogValue1);
    }

    // Check if there's serial input
    if (Serial.available() > 0) {
      String input = Serial.readStringUntil('\n');
      input.trim();  // Remove whitespace and newline characters

      if (input.equalsIgnoreCase("done")) {
        keepReading = false;
        Serial.println("Stopped reading analog values.");
      }
      // automatically set appropriate gamma correction
      if (analogValue0 < 420) {
        currentChALUT = ChA2LUT;
        currentChBLUT = ChA2LUT;
        Serial.print(F("LUT 2 SELECTED"));
        Serial.print("\n");
      } else {
        currentChALUT = ChA1LUT;
        currentChBLUT = ChA1LUT;
        Serial.print(F("LUT 1 SELECTED"));
        Serial.print("\n");
      }
    }
  }
}


///////////////////////////////////////// BIT REGISTERS ///////////////////////////////////////

////////////////////// TIMER 1 PWM FREQUENCY CONTROL //////////////////////////////

// Function to calculate the required prescaler and TOP value
long calculatePrescalerAndTOP(long desiredFrequency, long &prescaler) {
  long TOP = 0;

  // Possible prescaler values: 1, 8, 64, 256, 1024
  long possiblePrescalers[] = { 1, 8, 64, 256, 1024 };

  // Try each prescaler and calculate the corresponding TOP
  for (int i = 0; i < 5; i++) {
    long currentPrescaler = possiblePrescalers[i];

    // Calculate the TOP value
    long calculatedTOP = (CLOCK_FREQ / (2 * currentPrescaler * desiredFrequency)) - 1;

    // Check if the calculated TOP value is within the 16-bit range (0 to 65535)
    if (calculatedTOP >= 0 && calculatedTOP <= 65535) {
      prescaler = currentPrescaler;
      TOP = calculatedTOP;
      break;  // Stop after finding the first valid prescaler and TOP
    }
  }

  return TOP;
}

// Function to configure Timer1 with the calculated prescaler and TOP value
void configureTimer1(long prescaler, long TOP) {
  TCCR1A = 0;
  TCCR1B = 0;

  // Set Timer1 in 16-bit Phase Correct PWM mode
  TCCR1A |= (1 << COM1A1);  // Enable PWM on pin 9 (Channel A)
  TCCR1A |= (1 << COM1B1);  // Enable PWM on pin 10 (Channel B)
  TCCR1B |= (1 << WGM13);   // Set WGM13 bit
  TCCR1B &= ~(1 << WGM12);  // Clear WGM12 bit
  TCCR1A &= ~(1 << WGM11);  // Clear WGM11 bit
  TCCR1A &= ~(1 << WGM10);  // Clear WGM10 bit

  // Set the prescaler
  switch (prescaler) {
    case 1:
      TCCR1B |= (1 << CS10);  // Prescaler = 1
      break;
    case 8:
      TCCR1B |= (1 << CS11);  // Prescaler = 8
      break;
    case 64:
      TCCR1B |= (1 << CS11) | (1 << CS10);  // Prescaler = 64
      break;
    case 256:
      TCCR1B |= (1 << CS12);  // Prescaler = 256
      break;
    case 1024:
      TCCR1B |= (1 << CS12) | (1 << CS10);  // Prescaler = 1024
      break;
    default:
      break;
  }

  // Set the TOP value
  ICR1 = TOP;
}

// Function to constrain the frequency to a multiple of TABLE_SIZE
long constrainFrequency(long frequency) {
  return frequency - (frequency % TABLE_SIZE);
}


////////////////////////////// Timer3 interrupt control /////////////////////////////////

// Define a function pointer for the interrupt handler
void (*timer3Callback)() = nullptr;  // Initialize to null

// Function to set the callback for Timer3
void setTimer3Callback(void (*callback)()) {
  timer3Callback = callback;
}

// ISR for Timer3 Compare Match A - run the function each interrupt
ISR(TIMER3_COMPA_vect) {
  if (timer3Callback) {
    timer3Callback();  // Call the assigned callback function
  }
}

// Set timer3 frequency and configure for interrupts using CTC mode
void configureTimer3Interrupt(float frequency) {

  long t = micros();
  long prescaler = 0;
  long compareValue = 0;

  // Prescaler options: 1, 8, 64, 256, 1024
  long prescalerOptions[] = { 1, 8, 64, 256, 1024 };
  int prescalerIndex = 0;

  // Iterate through prescaler options to find a valid one
  for (prescalerIndex = 0; prescalerIndex < 5; prescalerIndex++) {
    prescaler = prescalerOptions[prescalerIndex];
    compareValue = (CLOCK_FREQ / (prescaler * frequency)) - 1;

    // Check if compareValue is within the valid 16-bit range
    if (compareValue >= 0 && compareValue <= 65535) {
      break;  // Found a valid prescaler and compareValue
    }
  }

  // If no valid prescaler is found, set to maximum possible values
  if (compareValue < 0 || compareValue > 65535) {
    Serial.println("Unable to configure timer for requested frequency. Adjusting to closest possible.");
    prescaler = 1024;
    compareValue = 65535;
  }

  if (compareValue < 0 || compareValue > 65535) {
    Serial.println("compare value bug!");
    compareValue = 65535;  // Ensure it fits in 16 bits
  }

  // Set Timer3 to CTC mode
  TCCR3A = 0;             // Normal operation
  TCCR3B = (1 << WGM32);  // CTC mode (clear on compare match)

  // Set the prescaler
  TCCR3B &= ~(1 << CS32 | 1 << CS31 | 1 << CS30);  // Clear prescaler bits
  switch (prescaler) {
    case 1: TCCR3B |= (1 << CS30); break;
    case 8: TCCR3B |= (1 << CS31); break;
    case 64: TCCR3B |= (1 << CS31) | (1 << CS30); break;
    case 256: TCCR3B |= (1 << CS32); break;
    case 1024: TCCR3B |= (1 << CS32) | (1 << CS30); break;
  }

  // reset the timer counter
  TCNT3 = 0;
  // Set the compare match value
  OCR3A = compareValue;

  // Clear the Timer3 Compare Match A interrupt flag
  TIFR3 |= (1 << OCF3A);  // Writing a 1 clears the flag

  // Enable Timer3 Compare Match A interrupt
  TIMSK3 |= (1 << OCIE3A);

  // Enable global interrupts
  sei();

  //t = micros() - t;
  //Serial.print("micros: ");
  //Serial.println(t);
}

// stop the interrupt function
void stopTimer3Interrupt() {
  // Disable the Timer3 interrupt
  TIMSK3 &= ~(1 << OCIE3A);  // Disable Timer3 Compare Match A interrupt
}
