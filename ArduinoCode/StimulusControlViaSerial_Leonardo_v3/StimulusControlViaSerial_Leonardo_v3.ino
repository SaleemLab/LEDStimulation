#define CLOCK_FREQ 16000000  // Arduino Leonardo clock frequency (16 MHz)
#define TABLE_SIZE 256       // Number of samples in the wavetable

// gamma correction LUTs

// Declare a pointer to the current LUT
const uint16_t* currentLUT;
char GammaLUTName;

const uint16_t PROGMEM defaultLUT[1041] = {
 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
 20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
 30, 31, 32, 33, 34, 35, 36, 37, 38, 39,
 40, 41, 42, 43, 44, 45, 46, 47, 48, 49,
 50, 51, 52, 53, 54, 55, 56, 57, 58, 59,
 60, 61, 62, 63, 64, 65, 66, 67, 68, 69,
 70, 71, 72, 73, 74, 75, 76, 77, 78, 79,
 80, 81, 82, 83, 84, 85, 86, 87, 88, 89,
 90, 91, 92, 93, 94, 95, 96, 97, 98, 99,
 100, 101, 102, 103, 104, 105, 106, 107, 108, 109,
 110, 111, 112, 113, 114, 115, 116, 117, 118, 119,
 120, 121, 122, 123, 124, 125, 126, 127, 128, 129,
 130, 131, 132, 133, 134, 135, 136, 137, 138, 139,
 140, 141, 142, 143, 144, 145, 146, 147, 148, 149,
 150, 151, 152, 153, 154, 155, 156, 157, 158, 159,
 160, 161, 162, 163, 164, 165, 166, 167, 168, 169,
 170, 171, 172, 173, 174, 175, 176, 177, 178, 179,
 180, 181, 182, 183, 184, 185, 186, 187, 188, 189,
 190, 191, 192, 193, 194, 195, 196, 197, 198, 199,
 200, 201, 202, 203, 204, 205, 206, 207, 208, 209,
 210, 211, 212, 213, 214, 215, 216, 217, 218, 219,
 220, 221, 222, 223, 224, 225, 226, 227, 228, 229,
 230, 231, 232, 233, 234, 235, 236, 237, 238, 239,
 240, 241, 242, 243, 244, 245, 246, 247, 248, 249,
 250, 251, 252, 253, 254, 255, 256, 257, 258, 259,
 260, 261, 262, 263, 264, 265, 266, 267, 268, 269,
 270, 271, 272, 273, 274, 275, 276, 277, 278, 279,
 280, 281, 282, 283, 284, 285, 286, 287, 288, 289,
 290, 291, 292, 293, 294, 295, 296, 297, 298, 299,
 300, 301, 302, 303, 304, 305, 306, 307, 308, 309,
 310, 311, 312, 313, 314, 315, 316, 317, 318, 319,
 320, 321, 322, 323, 324, 325, 326, 327, 328, 329,
 330, 331, 332, 333, 334, 335, 336, 337, 338, 339,
 340, 341, 342, 343, 344, 345, 346, 347, 348, 349,
 350, 351, 352, 353, 354, 355, 356, 357, 358, 359,
 360, 361, 362, 363, 364, 365, 366, 367, 368, 369,
 370, 371, 372, 373, 374, 375, 376, 377, 378, 379,
 380, 381, 382, 383, 384, 385, 386, 387, 388, 389,
 390, 391, 392, 393, 394, 395, 396, 397, 398, 399,
 400, 401, 402, 403, 404, 405, 406, 407, 408, 409,
 410, 411, 412, 413, 414, 415, 416, 417, 418, 419,
 420, 421, 422, 423, 424, 425, 426, 427, 428, 429,
 430, 431, 432, 433, 434, 435, 436, 437, 438, 439,
 440, 441, 442, 443, 444, 445, 446, 447, 448, 449,
 450, 451, 452, 453, 454, 455, 456, 457, 458, 459,
 460, 461, 462, 463, 464, 465, 466, 467, 468, 469,
 470, 471, 472, 473, 474, 475, 476, 477, 478, 479,
 480, 481, 482, 483, 484, 485, 486, 487, 488, 489,
 490, 491, 492, 493, 494, 495, 496, 497, 498, 499,
 500, 501, 502, 503, 504, 505, 506, 507, 508, 509,
 510, 511, 512, 513, 514, 515, 516, 517, 518, 519,
 520, 521, 522, 523, 524, 525, 526, 527, 528, 529,
 530, 531, 532, 533, 534, 535, 536, 537, 538, 539,
 540, 541, 542, 543, 544, 545, 546, 547, 548, 549,
 550, 551, 552, 553, 554, 555, 556, 557, 558, 559,
 560, 561, 562, 563, 564, 565, 566, 567, 568, 569,
 570, 571, 572, 573, 574, 575, 576, 577, 578, 579,
 580, 581, 582, 583, 584, 585, 586, 587, 588, 589,
 590, 591, 592, 593, 594, 595, 596, 597, 598, 599,
 600, 601, 602, 603, 604, 605, 606, 607, 608, 609,
 610, 611, 612, 613, 614, 615, 616, 617, 618, 619,
 620, 621, 622, 623, 624, 625, 626, 627, 628, 629,
 630, 631, 632, 633, 634, 635, 636, 637, 638, 639,
 640, 641, 642, 643, 644, 645, 646, 647, 648, 649,
 650, 651, 652, 653, 654, 655, 656, 657, 658, 659,
 660, 661, 662, 663, 664, 665, 666, 667, 668, 669,
 670, 671, 672, 673, 674, 675, 676, 677, 678, 679,
 680, 681, 682, 683, 684, 685, 686, 687, 688, 689,
 690, 691, 692, 693, 694, 695, 696, 697, 698, 699,
 700, 701, 702, 703, 704, 705, 706, 707, 708, 709,
 710, 711, 712, 713, 714, 715, 716, 717, 718, 719,
 720, 721, 722, 723, 724, 725, 726, 727, 728, 729,
 730, 731, 732, 733, 734, 735, 736, 737, 738, 739,
 740, 741, 742, 743, 744, 745, 746, 747, 748, 749,
 750, 751, 752, 753, 754, 755, 756, 757, 758, 759,
 760, 761, 762, 763, 764, 765, 766, 767, 768, 769,
 770, 771, 772, 773, 774, 775, 776, 777, 778, 779,
 780, 781, 782, 783, 784, 785, 786, 787, 788, 789,
 790, 791, 792, 793, 794, 795, 796, 797, 798, 799,
 800, 801, 802, 803, 804, 805, 806, 807, 808, 809,
 810, 811, 812, 813, 814, 815, 816, 817, 818, 819,
 820, 821, 822, 823, 824, 825, 826, 827, 828, 829,
 830, 831, 832, 833, 834, 835, 836, 837, 838, 839,
 840, 841, 842, 843, 844, 845, 846, 847, 848, 849,
 850, 851, 852, 853, 854, 855, 856, 857, 858, 859,
 860, 861, 862, 863, 864, 865, 866, 867, 868, 869,
 870, 871, 872, 873, 874, 875, 876, 877, 878, 879,
 880, 881, 882, 883, 884, 885, 886, 887, 888, 889,
 890, 891, 892, 893, 894, 895, 896, 897, 898, 899,
 900, 901, 902, 903, 904, 905, 906, 907, 908, 909,
 910, 911, 912, 913, 914, 915, 916, 917, 918, 919,
 920, 921, 922, 923, 924, 925, 926, 927, 928, 929,
 930, 931, 932, 933, 934, 935, 936, 937, 938, 939,
 940, 941, 942, 943, 944, 945, 946, 947, 948, 949,
 950, 951, 952, 953, 954, 955, 956, 957, 958, 959,
 960, 961, 962, 963, 964, 965, 966, 967, 968, 969,
 970, 971, 972, 973, 974, 975, 976, 977, 978, 979,
 980, 981, 982, 983, 984, 985, 986, 987, 988, 989,
 990, 991, 992, 993, 994, 995, 996, 997, 998, 999,
 1000, 1001, 1002, 1003, 1004, 1005, 1006, 1007, 1008, 1009,
 1010, 1011, 1012, 1013, 1014, 1015, 1016, 1017, 1018, 1019,
 1020, 1021, 1022, 1023, 1024, 1025, 1026, 1027, 1028, 1029,
 1030, 1031, 1032, 1033, 1034, 1035, 1036, 1037, 1038, 1039,
 1040
};

const uint16_t PROGMEM GREENLUT[1041] = {
    0, 3, 5, 8, 11, 14, 16, 19, 21, 23, 
    24, 25, 27, 28, 29, 30, 32, 33, 34, 36, 
    37, 38, 40, 41, 42, 43, 44, 45, 47, 48, 
    49, 50, 51, 52, 53, 54, 55, 56, 58, 59, 
    60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 
    70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 
    80, 82, 83, 84, 85, 86, 87, 88, 89, 90, 
    91, 92, 93, 94, 95, 96, 97, 99, 100, 101, 
    102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 
    112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 
    122, 122, 123, 124, 125, 126, 127, 128, 129, 130, 
    131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 
    141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 
    151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 
    161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 
    171, 172, 173, 174, 175, 176, 177, 177, 178, 179, 
    180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 
    190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 
    200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 
    210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 
    220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 
    230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 
    240, 241, 242, 243, 244, 245, 247, 248, 249, 250, 
    250, 251, 252, 253, 254, 255, 256, 257, 258, 259, 
    260, 260, 261, 262, 263, 264, 265, 266, 267, 268, 
    269, 270, 270, 271, 272, 273, 274, 275, 276, 277, 
    278, 279, 280, 281, 282, 283, 284, 285, 286, 287, 
    288, 289, 290, 291, 292, 293, 294, 295, 296, 297, 
    298, 299, 300, 301, 302, 302, 303, 304, 305, 306, 
    307, 308, 309, 310, 311, 312, 313, 314, 315, 316, 
    317, 318, 319, 319, 320, 321, 322, 323, 324, 325, 
    326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 
    335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 
    345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 
    355, 356, 357, 358, 359, 360, 360, 361, 362, 363, 
    364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 
    374, 375, 375, 376, 377, 378, 379, 380, 381, 382, 
    383, 384, 385, 386, 387, 388, 389, 389, 390, 391, 
    392, 393, 394, 395, 396, 397, 398, 399, 400, 401, 
    402, 403, 404, 405, 406, 407, 408, 408, 409, 410, 
    411, 412, 413, 414, 415, 416, 417, 418, 419, 420, 
    421, 422, 423, 424, 425, 426, 427, 428, 429, 430, 
    431, 432, 433, 434, 435, 436, 437, 438, 439, 440, 
    441, 442, 443, 443, 444, 445, 446, 447, 448, 449, 
    450, 451, 452, 453, 454, 455, 456, 456, 457, 458, 
    459, 460, 461, 462, 463, 464, 465, 466, 467, 468, 
    469, 470, 471, 472, 473, 474, 475, 476, 477, 478, 
    479, 480, 481, 482, 483, 484, 485, 486, 487, 488, 
    489, 490, 491, 492, 493, 494, 495, 496, 497, 498, 
    499, 500, 501, 502, 503, 504, 505, 506, 507, 508, 
    509, 510, 511, 512, 513, 514, 515, 516, 517, 518, 
    519, 520, 520, 521, 522, 523, 524, 525, 525, 526, 
    527, 528, 529, 530, 530, 531, 532, 533, 534, 535, 
    535, 536, 537, 538, 539, 540, 540, 541, 542, 543, 
    544, 545, 546, 547, 548, 549, 550, 551, 552, 553, 
    554, 555, 556, 556, 557, 558, 559, 560, 561, 562, 
    563, 564, 565, 566, 566, 567, 568, 569, 570, 571, 
    572, 573, 573, 574, 575, 576, 577, 578, 579, 579, 
    580, 581, 582, 583, 584, 585, 585, 586, 587, 588, 
    589, 590, 591, 591, 592, 593, 594, 595, 596, 596, 
    597, 598, 599, 600, 601, 602, 602, 603, 604, 605, 
    606, 607, 608, 609, 610, 610, 611, 612, 613, 614, 
    615, 616, 617, 618, 619, 620, 620, 621, 622, 623, 
    624, 625, 626, 627, 627, 628, 629, 630, 631, 632, 
    632, 633, 634, 635, 636, 637, 637, 638, 639, 640, 
    641, 642, 642, 643, 644, 645, 646, 647, 648, 649, 
    650, 651, 652, 653, 654, 655, 656, 657, 658, 659, 
    660, 661, 662, 663, 664, 665, 666, 667, 668, 668, 
    669, 670, 671, 672, 672, 673, 674, 675, 676, 676, 
    677, 678, 679, 680, 680, 681, 682, 683, 684, 684, 
    685, 686, 687, 688, 689, 690, 691, 692, 692, 693, 
    694, 695, 696, 697, 698, 699, 700, 701, 702, 703, 
    704, 705, 706, 707, 708, 709, 710, 711, 712, 713, 
    714, 715, 716, 717, 718, 719, 720, 721, 722, 723, 
    724, 725, 725, 726, 727, 728, 729, 730, 731, 731, 
    732, 733, 734, 734, 735, 736, 737, 737, 738, 739, 
    740, 740, 741, 742, 743, 743, 744, 745, 746, 746, 
    747, 748, 749, 749, 750, 751, 752, 753, 754, 755, 
    756, 756, 757, 758, 759, 760, 761, 762, 762, 763, 
    764, 765, 766, 767, 768, 768, 769, 770, 771, 772, 
    773, 774, 775, 776, 777, 777, 778, 779, 780, 781, 
    782, 783, 784, 785, 786, 787, 787, 788, 789, 790, 
    791, 792, 793, 794, 794, 795, 796, 797, 798, 799, 
    800, 800, 801, 802, 803, 804, 805, 806, 806, 807, 
    808, 809, 810, 811, 812, 813, 813, 814, 815, 816, 
    817, 818, 819, 820, 820, 821, 822, 823, 824, 825, 
    826, 827, 827, 828, 829, 830, 831, 832, 833, 834, 
    834, 835, 836, 837, 838, 839, 840, 840, 841, 842, 
    843, 844, 845, 846, 847, 847, 848, 849, 850, 851, 
    852, 853, 854, 855, 856, 857, 858, 859, 860, 861, 
    862, 863, 864, 865, 866, 867, 868, 869, 870, 871, 
    872, 873, 874, 874, 875, 876, 876, 877, 878, 879, 
    879, 880, 881, 881, 882, 883, 883, 884, 885, 886, 
    886, 887, 888, 888, 889, 890, 891, 891, 892, 893, 
    893, 894, 895, 896, 897, 898, 899, 900, 901, 901, 
    902, 903, 904, 905, 906, 907, 908, 909, 910, 911, 
    912, 913, 914, 915, 916, 916, 917, 918, 919, 920, 
    921, 922, 923, 924, 925, 926, 927, 928, 929, 930, 
    930, 931, 932, 933, 934, 935, 936, 937, 938, 939, 
    940, 941, 942, 942, 943, 944, 945, 946, 947, 948, 
    949, 950, 951, 952, 953, 953, 954, 955, 956, 957, 
    958, 959, 961, 962, 963, 964, 965, 966, 967, 968, 
    969, 971, 972, 973, 974, 975, 976, 977, 979, 980, 
    982, 984, 985, 987, 988, 990, 992, 993, 995, 997, 
    998};


const uint16_t PROGMEM UVLUT[1041] = {
    0, 0, 0, 0, 0, 0, 0, 0, 1, 3, // 1, 3 for last 2
    4, 5, 7, 8, 10, 11, 12, 14, 15, 17, 
    18, 19, 21, 22, 23, 24, 25, 27, 28, 29, 
    30, 31, 32, 33, 35, 36, 37, 38, 39, 40, 
    42, 43, 44, 45, 46, 47, 48, 50, 51, 52, 
    53, 54, 55, 57, 58, 59, 60, 61, 62, 64, 
    65, 66, 67, 68, 69, 70, 71, 73, 74, 75, 
    76, 77, 78, 79, 81, 82, 83, 84, 85, 86, 
    87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 
    97, 98, 99, 101, 102, 103, 104, 105, 106, 107, 
    108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 
    118, 120, 121, 122, 123, 124, 125, 126, 127, 128, 
    129, 130, 131, 132, 133, 134, 135, 137, 138, 139, 
    140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 
    150, 151, 152, 153, 155, 156, 157, 158, 159, 160, 
    161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 
    172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 
    182, 183, 185, 186, 187, 188, 189, 190, 191, 192, 
    193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 
    203, 204, 205, 207, 208, 209, 210, 211, 212, 213, 
    214, 215, 216, 217, 219, 220, 221, 222, 223, 224, 
    225, 226, 227, 229, 230, 231, 232, 232, 233, 234, 
    235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 
    245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 
    255, 257, 258, 259, 260, 261, 262, 263, 264, 265, 
    266, 267, 268, 269, 270, 272, 273, 274, 275, 276, 
    277, 278, 279, 280, 281, 282, 283, 284, 285, 286, 
    287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 
    297, 298, 299, 300, 301, 302, 303, 304, 305, 306, 
    307, 308, 309, 310, 311, 312, 313, 314, 315, 316, 
    317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 
    327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 
    337, 338, 339, 340, 342, 343, 344, 345, 346, 347, 
    348, 349, 350, 351, 352, 353, 354, 355, 356, 357, 
    358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 
    368, 369, 370, 371, 372, 373, 374, 375, 376, 377, 
    378, 379, 380, 381, 382, 383, 384, 385, 386, 387, 
    388, 390, 391, 392, 393, 394, 395, 396, 397, 398, 
    399, 400, 401, 402, 403, 404, 405, 406, 407, 408, 
    409, 410, 411, 413, 414, 415, 416, 417, 418, 419, 
    420, 421, 421, 422, 423, 424, 425, 426, 427, 428, 
    429, 430, 431, 432, 433, 434, 435, 436, 437, 438, 
    439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 
    449, 450, 451, 452, 454, 455, 456, 457, 458, 459, 
    460, 461, 462, 463, 464, 465, 466, 467, 468, 469, 
    471, 472, 473, 474, 475, 476, 477, 478, 479, 480, 
    481, 483, 484, 485, 486, 487, 488, 489, 490, 492, 
    493, 494, 495, 496, 497, 498, 499, 500, 501, 502, 
    503, 504, 505, 506, 507, 508, 509, 509, 510, 511, 
    512, 513, 514, 515, 516, 517, 518, 519, 520, 521, 
    521, 522, 523, 524, 525, 526, 527, 528, 529, 530, 
    531, 532, 533, 534, 535, 536, 537, 538, 539, 539, 
    540, 541, 542, 543, 544, 545, 546, 547, 548, 549, 
    550, 550, 551, 552, 553, 554, 555, 556, 557, 558, 
    559, 560, 561, 561, 562, 563, 564, 565, 566, 567, 
    568, 569, 570, 571, 572, 573, 574, 574, 575, 576, 
    577, 578, 579, 580, 581, 582, 583, 584, 585, 586, 
    587, 588, 589, 590, 591, 592, 593, 594, 595, 595, 
    596, 597, 598, 599, 600, 601, 602, 603, 604, 605, 
    606, 607, 608, 609, 610, 611, 612, 613, 614, 615, 
    615, 616, 617, 618, 619, 620, 621, 622, 623, 624, 
    625, 626, 627, 628, 629, 630, 631, 632, 633, 634, 
    635, 636, 637, 638, 639, 640, 641, 642, 643, 644, 
    645, 646, 647, 648, 649, 650, 650, 651, 652, 653, 
    654, 655, 656, 657, 658, 659, 660, 660, 661, 662, 
    663, 664, 665, 666, 667, 668, 669, 670, 671, 672, 
    673, 674, 676, 677, 678, 679, 680, 681, 682, 683, 
    684, 685, 686, 687, 688, 689, 690, 691, 692, 693, 
    694, 695, 695, 696, 697, 698, 699, 700, 701, 702, 
    703, 704, 705, 706, 707, 707, 708, 709, 710, 711, 
    712, 713, 713, 714, 715, 716, 717, 718, 719, 719, 
    720, 721, 722, 723, 724, 724, 725, 726, 727, 728, 
    729, 730, 731, 731, 732, 733, 734, 735, 736, 737, 
    738, 739, 740, 741, 741, 742, 743, 744, 745, 746, 
    747, 748, 749, 750, 751, 752, 753, 754, 755, 756, 
    757, 757, 758, 759, 760, 761, 762, 763, 764, 765, 
    766, 767, 768, 769, 770, 771, 772, 773, 774, 775, 
    776, 778, 779, 780, 781, 782, 783, 784, 785, 786, 
    787, 788, 789, 790, 791, 792, 793, 794, 794, 795, 
    796, 797, 798, 799, 800, 801, 802, 802, 803, 804, 
    805, 806, 807, 808, 809, 810, 810, 811, 812, 813, 
    814, 815, 816, 817, 818, 819, 820, 821, 821, 822, 
    823, 824, 825, 826, 827, 828, 829, 830, 831, 832, 
    833, 834, 835, 836, 837, 838, 839, 840, 841, 842, 
    843, 844, 845, 846, 847, 848, 849, 850, 851, 852, 
    853, 853, 854, 855, 856, 857, 857, 858, 859, 860, 
    860, 861, 862, 863, 864, 864, 865, 866, 867, 867, 
    868, 869, 870, 871, 871, 872, 873, 874, 875, 876, 
    877, 878, 879, 879, 880, 881, 882, 883, 884, 885, 
    886, 887, 888, 889, 890, 891, 892, 893, 894, 895, 
    896, 897, 898, 899, 900, 901, 902, 903, 904, 905, 
    906, 907, 908, 909, 910, 912, 913, 914, 915, 916, 
    917, 918, 919, 920, 921, 922, 923, 924, 925, 926, 
    927, 928, 929, 930, 931, 932, 933, 934, 935, 936, 
    937, 938, 939, 940, 942, 943, 944, 945, 946, 947, 
    948, 949, 950, 951, 952, 953, 954, 955, 956, 957, 
    959, 961, 963, 965, 967, 968, 970, 972, 974, 976, 
    982, 999, 999, 1000, 1000, 1001, 1002, 1002, 1003, 1003, 
    1004, 1005, 1005, 1006, 1006, 1007, 1008, 1008, 1009, 1009, 
    1010, 1011, 1011, 1012, 1012, 1013, 1014, 1014, 1015, 1015, 
    1016, 1017, 1017, 1018, 1018, 1019, 1020, 1020, 1021, 1021, 
    1022, 1023, 1023, 1024, 1024, 1025, 1026, 1026, 1027, 1027, 
    1028, 1029, 1029, 1030, 1030, 1031, 1032, 1032, 1033, 1033, 
    1034, 1035, 1035, 1036, 1036, 1037, 1038, 1038, 1039, 1039, 
    1040};


long prescaler;
long TOP;      // set by the clock
long TopLumi;  // use to limit max luminance
long desiredPWMFrequency = 7680;
float dutyCycle;

// serial
const byte numChars = 30;
char receivedChars[numChars];  // an array to store the received data
bool newData = false;

// stimulus selection char
String FirstChar;

// Array to hold the wavetable
uint16_t sineWaveTable[TABLE_SIZE];
//uint16_t baseSineWaveTable[TABLE_SIZE];

// wavetable step size
volatile int stepSize = 1;
volatile int tableIndex = 0; //for sinewave table
volatile int tableEnvIndex = 0; // for sinewave table contrast envelope
volatile int envCount = 0;  // contrast envelope counter
// contrast-envelope counter and contrast multiplier
volatile int nEnvCounts = 1;
volatile float contrastMult = 1;

// white noise update time and PWM value
volatile long updateTime;
volatile long randNumber = 0;

// flicker state
volatile bool toggleState = false;  // Flag to track the current state

// array for frozen white noise values
volatile uint16_t frozenWhiteNoiseTable[375];  // max number of values
volatile int frozenWhiteNoiseTableSize;



void setup() {
  pinMode(9, OUTPUT);      // Pin 9 controlled by Timer1 (Channel A)
  pinMode(4, OUTPUT);      // Pin 4 indicator pin, e.g. for sinewave cycles
  pinMode(5, OUTPUT);      // Pin 5 stim ON or OFF pin 

  PORTD &= ~(1 << PIND4);  // Ensure Pin 4 is set to LOW by changing register directly
  PORTC &= ~(1 << PORTC6); // Ensure Pin 5 is set to LOW


  Serial.begin(115200);

  // Set the desired PWM frequency (in Hz) -> set outside now
  //long desiredPWMFrequency = 7680;  // Example: 7680 Hz

  // Constrain the desired PWM frequency to be a multiple of TABLE_SIZE
  // (probably not important since precision seems low at high frequenies)
  desiredPWMFrequency = constrainFrequency(desiredPWMFrequency);

  // Calculate the prescaler and TOP value for the constrained frequency
  TOP = calculatePrescalerAndTOP(desiredPWMFrequency, prescaler);

  // Apply the prescaler and TOP value to Timer1
  configureTimer1(prescaler, TOP);

  TopLumi = TOP;  // set default TopLumi (max duty cycle) as the actual TOP (i.e. 100% for now)
  

  currentLUT = defaultLUT;
  GammaLUTName = 'd';

  // Generate the sine wave LUT based on the Timer1 config and TopLumi
  generateSineWaveTable(TopLumi);

  // initialise random number to 50% duty cyle for white noise stimuli
  randNumber = TopLumi / 2;

  // Set pin 9 to 50% duty cycle as default
  OCR1A = TopLumi / 2;

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
  Serial.println(receivedChars);
  char delimiters[] = ",";
  char *token;
  uint8_t idx = 0;
#define MAX_VALS 5  // max required? freq, duration, contrast, carrier freq?
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
    //Serial.println("Stim: Sinusoidal dimming");
    //Serial.flush();
    //Serial.print("Stim duration: ");
    //Serial.println(stimulusDuration);
    //Serial.flush();
    //Serial.print("Frequency: ");
    //Serial.println(frequency);
    //Serial.flush();
    
    outputSinewave(frequency, stimulusDuration);
  } else if (FirstChar == "wn")  // white noise
  {
    long stimulusDuration = atof(serialVals[1]);
    long updateTime = atof(serialVals[2]);
    //Serial.println("Stim: White noise");
    //Serial.print("Stim duration: ");
   // Serial.println(stimulusDuration);
    //Serial.flush();
   // Serial.print("Update time: ");
  // Serial.println(updateTime);
  //  Serial.flush();

    whiteNoise(updateTime, stimulusDuration);
  } else if (FirstChar == "fwn")  // frozen white noise
  {
    long stimulusDuration = atof(serialVals[1]);
    int updateTime = atof(serialVals[2]);
    int nReps = atof(serialVals[3]);
    //Serial.println("Stim: White noise");
    //Serial.flush();
  //  Serial.print("Stim duration: ");
  //  Serial.println(stimulusDuration);
  //  Serial.flush();
  //  Serial.print("Update time: ");
 //   Serial.println(updateTime);
  //  Serial.flush();

    frozenWhiteNoise(updateTime, stimulusDuration,nReps);
  } else if (FirstChar == "se")  // sinusoidal flicker with contrast envelope
  {
    long stimulusDuration = atof(serialVals[1]);
    float frequency = atof(serialVals[2]);
    float envFrequency = atof(serialVals[3]);
    //Serial.println("Stim: Sinusoidal env");
    //Serial.flush();
 //   Serial.print("Stim duration: ");
//    Serial.println(stimulusDuration);
 //   Serial.flush();
 //   Serial.print("Frequency: ");
 //   Serial.println(frequency);
 //   Serial.flush();
 //   Serial.print("Envelope freq: ");
 //   Serial.println(envFrequency);
 //   Serial.flush();

    SineContrastConv(stimulusDuration, frequency, envFrequency);
  } else if (FirstChar == "f")  // Square wave flicker
  {
    long stimulusDuration = atof(serialVals[1]);
    float FlickerFreq = atof(serialVals[2]);
//    Serial.println("Stim: Square Wave flicker");
//    Serial.flush();
 //   Serial.print("Stim duration: ");
 //   Serial.println(stimulusDuration);
 //   Serial.flush();
 //   Serial.print("Frequency: ");
 //   Serial.println(FlickerFreq);
 //   Serial.flush();

    FlickerLED(FlickerFreq, stimulusDuration);
  } else if (FirstChar == "st")  // Set TopLuminance
  {
    float TopMultiplier = atof(serialVals[1]);

    SetTopLumi(TopMultiplier);

  } else if (FirstChar == "sd")  // Set duty cycle
  {
    dutyCycle = atof(serialVals[1]);

    setDutyCycle(dutyCycle, TopLumi);

  } else if (FirstChar == "gc") // do gamma correction routine
  {
    float stepSize = atof(serialVals[1]);
    float waitTime = atof(serialVals[2]);
    int nReps = atoi(serialVals[3]);
    cycleDutyCycles(stepSize, waitTime, nReps, TopLumi);
  } else if (FirstChar == "agc") // apply gamma correction
  { 
    int GammaToUse = atoi(serialVals[1]);
    applyGammaCorrection(GammaToUse);
    
  } else if (FirstChar == "stat")
  {
    getStatus();
  }
  else  // not valid stimulus code
  {
    Serial.print(FirstChar);
    Serial.println(" is an invalid stimulus code - make sure you are using carriage return line ending");
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
    float angle = (2.0 * PI * i) / TABLE_SIZE;                     // Angle in radians
    sineWaveTable[i] = (uint16_t)((sin(angle) + 1.0) * (TOP / 2.0));  // Scale to 0-TOP, un-corrected sinewave
    // use this value to index into LUT
  }
     // for (int i = 0; i < 256; i++) {
    //    sineWaveTable[i] = baseSineWaveTable[i]; //currentLUT[baseSineWaveTable[i]]; // LUT-corrected 
   // }
}

void outputSinewave(float sinewaveFrequency, long duration) {

  tableIndex = 0;  // Start at the beginning of the sine wave table

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

  // Configure timer3 interrupt to updateFrequency
  configureTimer3Interrupt(updateFrequency);


  long startTime = millis();  // Record the start time
  PORTC |= (1 << PORTC6);
  // set timer3 interrupt callback function to play the sinewave
  setTimer3Callback(sinewaveInterrupt);

  // Loop until the specified duration has elapsed
  while (millis() - startTime < duration) {
    delayMicroseconds(1);  //wait for time to end
  }

  stopTimer3Interrupt();  // finish playing sinewave
  PORTD &= ~(1 << PIND4);  // Ensure Pin 4 is set to LOW by changing register directly
  PORTC &= ~(1 << PORTC6); // Ensure Pin 5 is set to LOW
  Serial.println("-1");

  OCR1A = TopLumi / 2;
}

// sinewave interrupt function
void sinewaveInterrupt() {
  //static int tableIndex = 0;  // Start at the beginning of the sine wave table
  // Update PWM duty cycle with the next sine wave value
  OCR1A = sineWaveTable[tableIndex];
  //Serial.print(OCR1A);
  //Serial.print(',');

 // tableIndex = (tableIndex + stepSize) % TABLE_SIZE;
  //if (tableIndex == 0) {
 //   PORTD ^= (1 << PIND4);  // Toggle Pin 4 if tableIndex is 0
 // }

  // Update the table index (wrap around if necessary)
  tableIndex = tableIndex + stepSize;
  if (tableIndex >= TABLE_SIZE) {
    PORTD ^= (1 << PIND4);  // Toggle Pin 4 if sine wave cycle finished
  }
  tableIndex = tableIndex % TABLE_SIZE;
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

  tableIndex = 0; // start at beginning of sinewave table
  tableEnvIndex = 0; // start at beginning of sinewave table
  envCount = 0;  // contrast envelope counter


  // Recalculate the effective update interval based on the step size
  float updateInterval = baseUpdateInterval * stepSize;
  float updateFrequency = 1e6 / updateInterval;  // update frequency for timer3 interrupt

  // Configure timer3 interrupt to updateFrequency
  configureTimer3Interrupt(updateFrequency);


  // now get interval for updating envelope
  float baseEnvUpdateInterval = 1.0 / (envelopeFreq * TABLE_SIZE);  // Time per table update in seconds
  baseEnvUpdateInterval *= 1e6;                                     // Convert to microseconds

  nEnvCounts = baseEnvUpdateInterval / updateInterval;
  //Serial.print("nc: ");
  //Serial.println(nEnvCounts);

  long startTime = millis();  // Record the start time
  // set timer3 interrupt callback function to play the sinewave
  PORTC |= (1 << PORTC6);  // Stim on pin 5
  setTimer3Callback(sinewaveEnvelopeInterrupt);

  // Loop until the specified duration has elapsed
  while (millis() - startTime < duration) {
    delayMicroseconds(1);  //wait for time to end
  }
  stopTimer3Interrupt();  // finish playing sinewave
  PORTD &= ~(1 << PIND4);  // Ensure Pin 4 is set to LOW by changing register directly
  PORTC &= ~(1 << PORTC6); // Ensure Pin 5 is set to LOW
  Serial.println("-1");
  Serial.flush();
  OCR1A = TopLumi / 2;
}



// sinewave contrast envelope interrupt function
void sinewaveEnvelopeInterrupt() {
  //static int tableIndex = 0;  // Start at the beginning of the sine wave table
  //static int tableEnvIndex = 0;
  //static int envCount = 0;  // contrast envelope counter

  // Update PWM duty cycle with the next sine wave value
  OCR1A = TopLumi / 2 + (sineWaveTable[tableIndex] - TopLumi / 2) * contrastMult; 
  
  // Update the table index (wrap around if necessary)
  //tableIndex = (tableIndex + stepSize) % TABLE_SIZE;
  //if (tableIndex == 0) {
  //  PORTD ^= (1 << PIND4);  // Toggle Pin 4 if tableIndex is 0
 // }

   tableIndex = tableIndex + stepSize;
  if (tableIndex >= TABLE_SIZE) {
    PORTD ^= (1 << PIND4);  // Toggle Pin 4 if sine wave cycle finished
  }
  tableIndex = tableIndex % TABLE_SIZE;

  // update counter for contrast envelope
  envCount = envCount + 1;
  if (envCount > nEnvCounts - 1)  // if time to update contrast envelope
  {
    envCount = 0;                                                  // reset to 1
    tableEnvIndex = (tableEnvIndex + 1) % TABLE_SIZE;              // get the next contrast value index
    contrastMult = sineWaveTable[tableEnvIndex] / float(TopLumi);  // use index to get the normalised contrast value
    //Serial.println(contrastMult);
  }
}


/////////////////////////////////// WHITE NOISE PWM FUNCTIONS //////////////////////////////

void whiteNoise(long updateTime, long duration) {

  float updateFrequency = 1e3 / updateTime;
  configureTimer3Interrupt(updateFrequency);


  Serial.print("TOP: ");
  Serial.println(TopLumi);

  long startTime = millis();  // Record the start time
  PORTC |= (1 << PORTC6);  // Stim on pin 5

  // set timer3 interrupt callback function to play the sinewave
  setTimer3Callback(whiteNoiseInterrupt);

  // Loop until the specified duration has elapsed
  while (millis() - startTime < duration) {
    delayMicroseconds(1);  //wait for time to end
  }

  stopTimer3Interrupt();  // stop white noise
  delay(updateTime);
  PORTD &= ~(1 << PIND4);  // Ensure Pin 4 is set to LOW by changing register directly
  PORTC &= ~(1 << PORTC6); // Ensure Pin 5 is set to LOW
  Serial.println("-1");
  Serial.flush();
  OCR1A = TopLumi / 2;
}

// whitenoise interrupt function
void whiteNoiseInterrupt() {
  // Update PWM duty cycle with the next sine wave value
  OCR1A = randNumber;
  PIND = (1 << PIND4);  // alternate PIN 4 value indicator pin
  Serial.println(randNumber);
  Serial.flush();
  randNumber = random(0, TopLumi);  // get a new random number ready
}


/////////////////////////////////// FROZEN WHITE NOISE PWM FUNCTIONS //////////////////////////////

void frozenWhiteNoise(int updateTime, long duration, long nReps) {

  float updateFrequency = 1e3 / updateTime;
  configureTimer3Interrupt(updateFrequency);

  frozenWhiteNoiseTableSize = duration / updateTime;
  long totalDuration = duration * nReps;
  Serial.print("LD: ");
  Serial.println(totalDuration);

  for (int i = 0; i < frozenWhiteNoiseTableSize; i++) {
    // Calculate the sine wave value (scaled between 0 and TOP)
    frozenWhiteNoiseTable[i] = (int)random(0, TopLumi);  // Scale to 0-TOP
  }

  Serial.print("TOP: ");
  Serial.println(TopLumi);
  
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
  PORTD &= ~(1 << PIND4);  // Ensure Pin 4 is set to LOW by changing register directly
  PORTC &= ~(1 << PORTC6); // Ensure Pin 5 is set to LOW
  Serial.println("-1");
  Serial.flush();
  OCR1A = TopLumi / 2;
}

// whitenoise interrupt function
void frozenWhiteNoiseInterrupt() {

  static int tableIndex = 0;  // Start at the beginning of the sine wave table
  // Update PWM duty cycle with the frozen white noise value
  OCR1A = frozenWhiteNoiseTable[tableIndex];
  PIND = (1 << PIND4);  // alternate PIN 4 value indicator pin

  //if (tableIndex == 0) {
  //  PORTD ^= (1 << PIND4);  // Toggle Pin 4 if tableIndex is 0
  //}
  //Serial.print("ti: ");
  //Serial.println(tableIndex);
  Serial.println(frozenWhiteNoiseTable[tableIndex]);
  Serial.flush();
  // Update the table index (wrap around at actual white noise table size)
  tableIndex = (tableIndex + 1) % frozenWhiteNoiseTableSize;  //
}

////////////////////////////////////////// SQUARE WAVE FLICKER STIMULUS ////////////////////////////////
void FlickerLED(float flickerFreq, long duration) {
  const long interval = 1e6 / (flickerFreq * 2);  // interval at which to blink (microseconds)
  float updateFrequency = 1e6 / interval;         // update frequency for timer3 interrupt

  configureTimer3Interrupt(updateFrequency);

  long startTime = millis();  // Record the start time
  PORTC |= (1 << PORTC6);  // Stim on pin 5
  setTimer3Callback(SquareWaveFlickerInterrupt);

  // Loop until the specified duration has elapsed
  while (millis() - startTime < duration) {
    delayMicroseconds(1);  //wait for time to end
  }

  stopTimer3Interrupt();  // stop flicker
  PORTD &= ~(1 << PIND4);  // Ensure Pin 4 is set to LOW by changing register directly
  PORTC &= ~(1 << PORTC6); // Ensure Pin 5 is set to LOW
  Serial.println("-1");
  Serial.flush();
  OCR1A = TopLumi / 2;
}

// square wave flicker interrupt function
void SquareWaveFlickerInterrupt() {
  OCR1A = (OCR1A == TopLumi) ? 0 : TopLumi;  // Toggle between 0 and TOP
  PORTD ^= (1 << PIND4);                     // set indicator pin similarly
}


/////////////////////////////////// SOME GENERIC PWM FUNCTIONS ///////////////////////////////////////////
// function to artifically lower the max PWM duty cycle. (i.e. TopMultiplier=0.5 means max duty cycle of 50%)
// other functions will work as normal but scale to this TOP value
void SetTopLumi(float TopMultiplier) {

  if (TopMultiplier > 1) {
    TopMultiplier = 1;
  } else if (TopMultiplier <= 0) {
    TopMultiplier = 1;
  }

  // get new TOP value to use
  TopLumi = float(TOP) * TopMultiplier;

  // Generate the sine wave LUT based on the Timer1 config
  generateSineWaveTable(TopLumi);

  // initialise random number to 50% duty cyle
  randNumber = TopLumi / 2;

  // Set pin 9 to 50% duty cycle as default
  OCR1A = TopLumi / 2;
}


// set the duty cycle manually until a new value is requested
void setDutyCycle(float dutyCyclePercentage, long TopLumi) {
  // Constrain the duty cycle percentage between 0% and 100%
  if (dutyCyclePercentage < 0.0) dutyCyclePercentage = 0.0;
  if (dutyCyclePercentage > 100.0) dutyCyclePercentage = 100.0;

  // Calculate the OCR1A value based on the duty cycle and TOP
  uint16_t ocrValue = (long)((dutyCyclePercentage / 100.0) * TopLumi);
  
  // Set OCR1A to control the duty cycle
  //OCR1A = ocrValue;
  OCR1A = pgm_read_word_near(currentLUT+ocrValue);
  Serial.print(ocrValue);
  Serial.print(',');
  Serial.println(OCR1A);
}


// Run through duty cycles to perform gamma correction
void cycleDutyCycles(float stepSize, float waitTime, int nReps, long TopLumi){
  float dutyCycle = 0;

  for (int irep=0; irep<nReps; irep++)
  {
  while (dutyCycle<=1)
  {
    Serial.println(dutyCycle);
    long ocrValue = (long)(dutyCycle* TopLumi);
    OCR1A = ocrValue;
    delay(waitTime);
    dutyCycle = dutyCycle+stepSize;
  }
  }
  Serial.println("-1");
  // Set pin 9 to 50% duty cycle as default
  OCR1A = TopLumi / 2;
    
}

void applyGammaCorrection(int LUT_index){
  if (LUT_index==0)
    {
      currentLUT = defaultLUT;
      // update sinewavetable!
      Serial.println(F("defaultLUT selected"));
      GammaLUTName = 'd';
    
    } 
    else if (LUT_index==1)
    {
      currentLUT = GREENLUT;
      Serial.println(F("GREENLUT selected"));
      GammaLUTName = 'g';

    }
    else if (LUT_index==2)
    {
      currentLUT = UVLUT;
      Serial.println(F("UVLUT selected"));
      GammaLUTName = 'u';

    }
    else
    {
      Serial.println("Invalid LUT_index selected, gammaLUT not changed.");
    }
}

void getStatus()
{
  Serial.print("PWM FREQ: ");
  Serial.println(desiredPWMFrequency);
  Serial.print("TOP: ");
  Serial.println(TOP);
  Serial.print("TopLumi: ");
  Serial.println(TopLumi);
  Serial.print("Duty cycle: ");
  Serial.println(dutyCycle);
  Serial.print("Gamma Correction LUT: ");
  Serial.println(GammaLUTName);

}


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
