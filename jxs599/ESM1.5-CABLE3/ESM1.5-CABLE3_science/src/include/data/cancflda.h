! CANCFLDA List of Ancillary Fields - Atmosphere Stash Codes,
! Model Codes and Logical file Numbers
!
! -------------------------------------------------------------------
      INTEGER :: ITEM_CODES_ANCIL(NANCIL_FIELDS)  ! Stash Codes
      INTEGER :: MODEL_CODES_ANCIL(NANCIL_FIELDS) ! Model Codes
      INTEGER :: ANCIL_FILE_NO(NANCIL_FIELDS)     ! Logical file numbers

! -----------------------------------------
! Note 127:154 not used in UM ; set to zero
! -----------------------------------------

      DATA ITEM_CODES_ANCIL(1:100)/                                     &
     &  30,  33,  34,  35,  36,  37,  60,   0,  23,  20,                &
     &  40,  41,   0,  43,  44,   0,  46,  47,  50,  51,                &
     &  52,  53,  54,  55,  56,  26,  31,  24,  32,  28,                &
     &  29,  93,   0,   0,  48,   9,   0,   0,  58,  59,                &
     &  88,  87,  85,  57,  90,  17,  18, 301, 302, 303,                &
     & 304, 305, 306, 307, 308, 309, 310, 311, 312, 313,                &
     & 314, 315, 316, 317, 318, 319, 320, 127, 128, 129,                &
     &   0, 121, 122, 123, 124, 125, 126, 251, 207, 208,                &
     & 209, 160, 216, 217, 218, 213, 219, 220, 223, 321,                &
     & 322, 323, 324, 325, 326, 327, 328, 329, 330, 331/

    !kdcorbin, 05/10 - added 20 tracer flux variables
      DATA ITEM_CODES_ANCIL(101:NANCIL_FIELDS)/                         &
     & 332, 333, 334, 335, 336, 337, 338, 339, 340, 341,                &
     & 505, 418, 419, 420, 421, 422, 423, 424, 425, 426,                &
     & 130, 131, 132, 153, 151, 152,   0,   0,   0,   0,                &
     &   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,                &
     &   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,                &
     &   0,   0,   0,   0,   5,   6, 351, 352, 353, 354,                &
     & 355, 356, 357, 358, 359, 360, 361, 362, 363, 364,                &
     & 365, 366, 367, 368, 369, 370, 371, 480, 481, 482,                &
     & 483, 484, 485, 486, 487, 134, 135,3100,3101,3102,                &
     &3103,3104,3105,3106,3107,3108,3109,3110,3111,3112,                &
     &3113,3114,3115,3116,3117,3118,3119,884/

      DATA MODEL_CODES_ANCIL(1:100) /                                   &
     &   1,   1,   1,   1,   1,   1,   1,   0,   1,   1,                &
     &   1,   1,   1,   0,   1,   0,   1,   1,   1,   1,                &
     &   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,                &
     &   1,   1,   0,   0,   1,   1,   0,   0,   1,   1,                &
     &   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,                &
     &   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,                &
     &   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,                &
     &   0,   1,   1,   1,   1,   1,   1,   1,   1,   1,                &
     &   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,                &
     &   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/

      DATA MODEL_CODES_ANCIL(101:NANCIL_FIELDS) /                       &
     &   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,                &
     &   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,                &
     &   1,   1,   1,   1,   1,   1,   0,   0,   0,   0,                &
     &   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,                &
     &   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,                &
     &   0,   0,   0,   0,   1,   1,   1,   1,   1,   1,                &
     &   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,                &
     &   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,                &
     &   1,   1,   1,   1,   1,   1,   1,   1/
 
      DATA ANCIL_FILE_NO(1:100) /                                       &
     &   9,  10,  10,  10,  10,  10,   1,   0,   2,   3,                &
     &   4,   4,   0,   4,   4,   0,   4,   4,   5,   5,                &
     &   5,   5,   5,   5,   5,   5,   7,   6,   7,   8,                &
     &   8,   9,   0,   0,   4,   2,   0,   0,  12,  12,                &
     &  13,  13,  13,  14,  14,  10,  10,  15,  15,  15,                &
     &  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,                &
     &  15,  15,  15,  15,  15,  15,  15,  12,  23,  23,                &
     &   0,  17,  18,  18,  18,  18,  12,  24,   4,   5,                &
     &   5,  19,  20,  21,  21,  21,  22,   4,   4,  16,                &
     &  16,  16,  16,  16,  16,  16,  16,  16,  16,  16/

    !kdcorbin, 05/10 - added 20 tracer flux variables
      DATA ANCIL_FILE_NO(101:NANCIL_FIELDS) /                           &
     &  16,  16,  16,  16,  16,  16,  16,  16,  16,  25,                &
     &  26,  27,  27,  27,  27,  27,  27,  27,  27,  27,                &
     &  28,  28,  29,  30,  31,  31,   0,   0,   0,   0,                &
     &   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,                &
     &   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,                &
     &   0,   0,   0,   0,  10,  10,  38,  39,  39,  39,                &
     &  40,  40,  41,  41,  42,  42,  42,  43,  43,  43,                &
     &  43,  43,  43,  44,  44,  44,  45,  46,  46,  46,                &
     &  46,  46,  46,  46,  46,  47,  47,  48,  49,  50,                &
     &  51,  52,  53,  54,  55,  56,  57,  58,  59,  60,                &
     &  61,  62,  63,  64,  65,  66,  67,  68/
! CANCFLDA end
