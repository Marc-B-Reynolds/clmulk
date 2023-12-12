// all primitive elements that fix in a 16-bit package
const uint16_t gf2_primitives_u16[] = {
  0x0007,0x000b,0x000d,0x0013,0x0019,0x0025,0x0029,0x002f,
  0x0037,0x003b,0x003d,0x0043,0x005b,0x0061,0x0067,0x006d,
  0x0073,0x0083,0x0089,0x008f,0x0091,0x009d,0x00a7,0x00ab,
  0x00b9,0x00bf,0x00c1,0x00cb,0x00d3,0x00d5,0x00e5,0x00ef,
  0x00f1,0x00f7,0x00fd,0x011d,0x012b,0x012d,0x014d,0x015f,
  0x0163,0x0165,0x0169,0x0171,0x0187,0x018d,0x01a9,0x01c3,
  0x01cf,0x01e7,0x01f5,0x0211,0x021b,0x0221,0x022d,0x0233,
  0x0259,0x025f,0x0269,0x026f,0x0277,0x027d,0x0287,0x0295,
  0x02a3,0x02a5,0x02af,0x02b7,0x02bd,0x02cf,0x02d1,0x02db,
  0x02f5,0x02f9,0x0313,0x0315,0x031f,0x0323,0x0331,0x033b,
  0x034f,0x035b,0x0361,0x036b,0x036d,0x0373,0x037f,0x0385,
  0x038f,0x03b5,0x03b9,0x03c7,0x03cb,0x03cd,0x03d5,0x03d9,
  0x03e3,0x03e9,0x03fb,0x0409,0x041b,0x0427,0x042d,0x0465,
  0x046f,0x0481,0x048b,0x04c5,0x04d7,0x04e7,0x04f3,0x04ff,
  0x050d,0x0519,0x0523,0x0531,0x053d,0x0543,0x0557,0x056b,
  0x0585,0x058f,0x0597,0x05a1,0x05c7,0x05e5,0x05f7,0x05fb,
  0x0613,0x0615,0x0625,0x0637,0x0643,0x064f,0x065b,0x0679,
  0x067f,0x0689,0x06b5,0x06c1,0x06d3,0x06df,0x06fd,0x0717,
  0x071d,0x0721,0x0739,0x0747,0x074d,0x0755,0x0759,0x0763,
  0x077d,0x078d,0x0793,0x07b1,0x07db,0x07f3,0x07f9,0x0805,
  0x0817,0x082b,0x082d,0x0847,0x0863,0x0865,0x0871,0x087b,
  0x088d,0x0895,0x089f,0x08a9,0x08b1,0x08cf,0x08d1,0x08e1,
  0x08e7,0x08eb,0x08f5,0x090d,0x0913,0x0925,0x0929,0x093b,
  0x093d,0x0945,0x0949,0x0951,0x095b,0x0973,0x0975,0x097f,
  0x0983,0x098f,0x09ab,0x09ad,0x09b9,0x09c7,0x09d9,0x09e5,
  0x09f7,0x0a01,0x0a07,0x0a13,0x0a15,0x0a29,0x0a49,0x0a61,
  0x0a6d,0x0a79,0x0a7f,0x0a85,0x0a91,0x0a9d,0x0aa7,0x0aab,
  0x0ab3,0x0ab5,0x0ad5,0x0adf,0x0ae9,0x0aef,0x0af1,0x0afb,
  0x0b03,0x0b09,0x0b11,0x0b33,0x0b3f,0x0b41,0x0b4b,0x0b59,
  0x0b5f,0x0b65,0x0b6f,0x0b7d,0x0b87,0x0b8b,0x0b93,0x0b95,
  0x0baf,0x0bb7,0x0bbd,0x0bc9,0x0bdb,0x0bdd,0x0be7,0x0bed,
  0x0c0b,0x0c0d,0x0c19,0x0c1f,0x0c57,0x0c61,0x0c6b,0x0c73,
  0x0c85,0x0c89,0x0c97,0x0c9b,0x0c9d,0x0cb3,0x0cbf,0x0cc7,
  0x0ccd,0x0cd3,0x0cd5,0x0ce3,0x0ce9,0x0cf7,0x0d03,0x0d0f,
  0x0d1d,0x0d27,0x0d2d,0x0d41,0x0d47,0x0d55,0x0d59,0x0d63,
  0x0d6f,0x0d71,0x0d93,0x0d9f,0x0da9,0x0dbb,0x0dbd,0x0dc9,
  0x0dd7,0x0ddb,0x0de1,0x0de7,0x0df5,0x0e05,0x0e1d,0x0e21,
  0x0e27,0x0e2b,0x0e33,0x0e39,0x0e47,0x0e4b,0x0e55,0x0e5f,
  0x0e71,0x0e7b,0x0e7d,0x0e81,0x0e93,0x0e9f,0x0ea3,0x0ebb,
  0x0ecf,0x0edd,0x0ef3,0x0ef9,0x0f0b,0x0f19,0x0f31,0x0f37,
  0x0f5d,0x0f6b,0x0f6d,0x0f75,0x0f83,0x0f91,0x0f97,0x0f9b,
  0x0fa7,0x0fad,0x0fb5,0x0fcd,0x0fd3,0x0fe5,0x0fe9,0x1053,
  0x1069,0x107b,0x107d,0x1099,0x10d1,0x10eb,0x1107,0x111f,
  0x1123,0x113b,0x114f,0x1157,0x1161,0x116b,0x1185,0x11b3,
  0x11d9,0x11df,0x120d,0x1237,0x123d,0x1267,0x1273,0x127f,
  0x12b9,0x12c1,0x12cb,0x130f,0x131d,0x1321,0x1339,0x133f,
  0x134d,0x1371,0x1399,0x13a3,0x13a9,0x1407,0x1431,0x1437,
  0x144f,0x145d,0x1467,0x1475,0x14a7,0x14ad,0x14d3,0x150f,
  0x151d,0x154d,0x1593,0x15c5,0x15d7,0x15dd,0x15eb,0x1609,
  0x1647,0x1655,0x1659,0x16a5,0x16bd,0x1715,0x1719,0x1743,
  0x1745,0x1775,0x1789,0x17ad,0x17b3,0x17bf,0x17c1,0x1857,
  0x185d,0x1891,0x1897,0x18b9,0x18ef,0x191b,0x1935,0x1941,
  0x1965,0x197b,0x198b,0x19b1,0x19bd,0x19c9,0x19cf,0x19e7,
  0x1a1b,0x1a2b,0x1a33,0x1a69,0x1a8b,0x1ad1,0x1ae1,0x1af5,
  0x1b0b,0x1b13,0x1b1f,0x1b57,0x1b91,0x1ba7,0x1bbf,0x1bc1,
  0x1bd3,0x1c05,0x1c11,0x1c17,0x1c27,0x1c4d,0x1c87,0x1c9f,
  0x1ca5,0x1cbb,0x1cc5,0x1cc9,0x1ccf,0x1cf3,0x1d07,0x1d23,
  0x1d43,0x1d51,0x1d5b,0x1d75,0x1d85,0x1d89,0x1e15,0x1e19,
  0x1e2f,0x1e45,0x1e51,0x1e67,0x1e73,0x1e8f,0x1ee3,0x1f11,
  0x1f1b,0x1f27,0x1f71,0x1f99,0x1fbb,0x1fbd,0x1fc9,0x201b,
  0x2027,0x2035,0x2053,0x2065,0x206f,0x208b,0x208d,0x209f,
  0x20a5,0x20af,0x20bb,0x20bd,0x20c3,0x20c9,0x20e1,0x20f3,
  0x210d,0x2115,0x2129,0x212f,0x213b,0x2143,0x2167,0x216b,
  0x2179,0x2189,0x2197,0x219d,0x21bf,0x21c1,0x21c7,0x21cd,
  0x21df,0x21e3,0x21f1,0x21fb,0x2219,0x2225,0x2237,0x223d,
  0x2243,0x225b,0x225d,0x2279,0x227f,0x2289,0x2297,0x229b,
  0x22b3,0x22bf,0x22cd,0x22ef,0x22f7,0x22fb,0x2305,0x2327,
  0x232b,0x2347,0x2355,0x2359,0x236f,0x2371,0x237d,0x2387,
  0x238d,0x2395,0x23a3,0x23a9,0x23b1,0x23b7,0x23bb,0x23e1,
  0x23ed,0x23f9,0x240b,0x2413,0x241f,0x2425,0x2429,0x243d,
  0x2451,0x2457,0x2461,0x246d,0x247f,0x2483,0x249b,0x249d,
  0x24b5,0x24bf,0x24c1,0x24c7,0x24cb,0x24e3,0x2509,0x2517,
  0x251d,0x2521,0x252d,0x2539,0x2553,0x2555,0x2563,0x2571,
  0x2577,0x2587,0x258b,0x2595,0x2599,0x259f,0x25af,0x25bd,
  0x25c5,0x25cf,0x25d7,0x25eb,0x2603,0x2605,0x2611,0x262d,
  0x263f,0x264b,0x2653,0x2659,0x2669,0x2677,0x267b,0x2687,
  0x2693,0x2699,0x26b1,0x26b7,0x26bd,0x26c3,0x26eb,0x26f5,
  0x2713,0x2729,0x273b,0x274f,0x2757,0x275d,0x276b,0x2773,
  0x2779,0x2783,0x2791,0x27a1,0x27b9,0x27c7,0x27cb,0x27df,
  0x27ef,0x27f1,0x2807,0x2819,0x281f,0x2823,0x2831,0x283b,
  0x283d,0x2845,0x2867,0x2875,0x2885,0x28ab,0x28ad,0x28bf,
  0x28cd,0x28d5,0x28df,0x28e3,0x28e9,0x28fb,0x2909,0x290f,
  0x2911,0x291b,0x292b,0x2935,0x293f,0x2941,0x294b,0x2955,
  0x2977,0x297d,0x2981,0x2993,0x299f,0x29af,0x29b7,0x29bd,
  0x29c3,0x29d7,0x29f3,0x29f5,0x2a03,0x2a0f,0x2a1d,0x2a21,
  0x2a33,0x2a35,0x2a4d,0x2a69,0x2a6f,0x2a71,0x2a7b,0x2a7d,
  0x2aa5,0x2aa9,0x2ab1,0x2ac5,0x2ad7,0x2adb,0x2aeb,0x2af3,
  0x2b01,0x2b15,0x2b23,0x2b25,0x2b2f,0x2b37,0x2b43,0x2b49,
  0x2b6d,0x2b7f,0x2b85,0x2b97,0x2b9b,0x2bad,0x2bb3,0x2bd9,
  0x2be5,0x2bfd,0x2c0f,0x2c21,0x2c2b,0x2c2d,0x2c3f,0x2c41,
  0x2c4d,0x2c71,0x2c8b,0x2c8d,0x2c95,0x2ca3,0x2caf,0x2cbd,
  0x2cc5,0x2cd1,0x2cd7,0x2ce1,0x2ce7,0x2ceb,0x2d0d,0x2d19,
  0x2d29,0x2d2f,0x2d37,0x2d3b,0x2d45,0x2d5b,0x2d67,0x2d75,
  0x2d89,0x2d8f,0x2da7,0x2dab,0x2db5,0x2de3,0x2df1,0x2dfd,
  0x2e07,0x2e13,0x2e15,0x2e29,0x2e49,0x2e4f,0x2e5b,0x2e5d,
  0x2e61,0x2e6b,0x2e8f,0x2e91,0x2e97,0x2e9d,0x2eab,0x2eb3,
  0x2eb9,0x2edf,0x2efb,0x2efd,0x2f05,0x2f09,0x2f11,0x2f17,
  0x2f3f,0x2f41,0x2f4b,0x2f4d,0x2f59,0x2f5f,0x2f65,0x2f69,
  0x2f95,0x2fa5,0x2faf,0x2fb1,0x2fcf,0x2fdd,0x2fe7,0x2fed,
  0x2ff5,0x2fff,0x3007,0x3015,0x3019,0x302f,0x3049,0x304f,
  0x3067,0x3079,0x307f,0x3091,0x30a1,0x30b5,0x30bf,0x30c1,
  0x30d3,0x30d9,0x30e5,0x30ef,0x3105,0x310f,0x3135,0x3147,
  0x314d,0x315f,0x3163,0x3171,0x317b,0x31a3,0x31a9,0x31b7,
  0x31c5,0x31c9,0x31db,0x31e1,0x31eb,0x31ed,0x31f3,0x31ff,
  0x3209,0x320f,0x321d,0x3227,0x3239,0x324b,0x3253,0x3259,
  0x3265,0x3281,0x3293,0x3299,0x329f,0x32a9,0x32b7,0x32bb,
  0x32c3,0x32d7,0x32db,0x32e7,0x3307,0x3315,0x332f,0x3351,
  0x335d,0x3375,0x3397,0x339b,0x33ab,0x33b9,0x33c1,0x33c7,
  0x33d5,0x33e3,0x33e5,0x33f7,0x33fb,0x3409,0x341b,0x3427,
  0x3441,0x344d,0x345f,0x3469,0x3477,0x347b,0x3487,0x3493,
  0x3499,0x34a5,0x34bd,0x34c9,0x34db,0x34e7,0x34f9,0x350d,
  0x351f,0x3525,0x3531,0x3537,0x3545,0x354f,0x355d,0x356d,
  0x3573,0x357f,0x359d,0x35a1,0x35b9,0x35cd,0x35d5,0x35d9,
  0x35e3,0x35e9,0x35ef,0x3601,0x360b,0x361f,0x3625,0x362f,
  0x363b,0x3649,0x3651,0x365b,0x3673,0x3675,0x3691,0x369b,
  0x369d,0x36ad,0x36cb,0x36d3,0x36d5,0x36e3,0x36ef,0x3705,
  0x370f,0x371b,0x3721,0x372d,0x3739,0x3741,0x3747,0x3753,
  0x3771,0x3777,0x378b,0x3795,0x3799,0x37a3,0x37c5,0x37cf,
  0x37d1,0x37d7,0x37dd,0x37e1,0x37f3,0x3803,0x3805,0x3817,
  0x381d,0x3827,0x3833,0x384b,0x3859,0x3869,0x3871,0x38a3,
  0x38b1,0x38bb,0x38c9,0x38cf,0x38e1,0x38f3,0x38f9,0x3901,
  0x3907,0x390b,0x3913,0x3931,0x394f,0x3967,0x396d,0x3983,
  0x3985,0x3997,0x39a1,0x39a7,0x39ad,0x39cb,0x39cd,0x39d3,
  0x39ef,0x39f7,0x39fd,0x3a07,0x3a29,0x3a2f,0x3a3d,0x3a51,
  0x3a5d,0x3a61,0x3a67,0x3a73,0x3a75,0x3a89,0x3ab9,0x3abf,
  0x3acd,0x3ad3,0x3ad5,0x3adf,0x3ae5,0x3ae9,0x3afb,0x3b11,
  0x3b2b,0x3b2d,0x3b35,0x3b3f,0x3b53,0x3b59,0x3b63,0x3b65,
  0x3b6f,0x3b71,0x3b77,0x3b8b,0x3b99,0x3ba5,0x3ba9,0x3bb7,
  0x3bbb,0x3bd1,0x3be7,0x3bf3,0x3bff,0x3c0d,0x3c13,0x3c15,
  0x3c1f,0x3c23,0x3c25,0x3c3b,0x3c4f,0x3c5d,0x3c6d,0x3c83,
  0x3c8f,0x3c9d,0x3ca7,0x3cab,0x3cb9,0x3cc7,0x3ce9,0x3cfb,
  0x3cfd,0x3d03,0x3d17,0x3d1b,0x3d21,0x3d2d,0x3d33,0x3d35,
  0x3d41,0x3d4d,0x3d65,0x3d69,0x3d7d,0x3d81,0x3d95,0x3db1,
  0x3db7,0x3dc3,0x3dd1,0x3ddb,0x3de7,0x3deb,0x3df9,0x3e05,
  0x3e09,0x3e0f,0x3e1b,0x3e2b,0x3e3f,0x3e41,0x3e53,0x3e65,
  0x3e69,0x3e8b,0x3ea3,0x3ebd,0x3ec5,0x3ed7,0x3edd,0x3ee1,
  0x3ef9,0x3f0d,0x3f19,0x3f1f,0x3f25,0x3f37,0x3f3d,0x3f43,
  0x3f45,0x3f49,0x3f51,0x3f57,0x3f61,0x3f83,0x3f89,0x3f91,
  0x3fab,0x3fb5,0x3fe3,0x3ff7,0x3ffd,0x402b,0x4039,0x4053,
  0x405f,0x407b,0x40a9,0x40af,0x40bb,0x40bd,0x40cf,0x40eb,
  0x40f3,0x410d,0x4113,0x413b,0x4143,0x419b,0x419d,0x41a7,
  0x41ad,0x41b5,0x41d5,0x41d9,0x41f1,0x420d,0x4257,0x4261,
  0x427f,0x4285,0x429d,0x42c7,0x42cb,0x42cd,0x42e3,0x42e9,
  0x42ef,0x4309,0x4321,0x433f,0x437d,0x4387,0x4395,0x43af,
  0x43c9,0x43eb,0x43ed,0x440b,0x4443,0x4473,0x44d3,0x44d5,
  0x44df,0x44e3,0x44fb,0x452b,0x4539,0x4559,0x456f,0x4599,
  0x459f,0x45a5,0x45b7,0x45c5,0x45d7,0x45e7,0x45f3,0x45ff,
  0x460f,0x461d,0x4627,0x4635,0x4647,0x4659,0x4663,0x4671,
  0x467b,0x46a5,0x46c5,0x46cf,0x46db,0x4731,0x474f,0x477f,
  0x47a7,0x47c1,0x47e5,0x47e9,0x47ef,0x4813,0x4819,0x483b,
  0x4843,0x485b,0x4861,0x4867,0x487f,0x4883,0x48e3,0x491d,
  0x492d,0x4939,0x4969,0x4977,0x4987,0x49a9,0x49cf,0x49e1,
  0x4a35,0x4a65,0x4a7b,0x4a81,0x4a87,0x4a8b,0x4aa3,0x4ac9,
  0x4add,0x4af9,0x4b37,0x4b49,0x4b5b,0x4b6b,0x4b6d,0x4b73,
  0x4b79,0x4b85,0x4b97,0x4ba1,0x4bab,0x4bb9,0x4bcb,0x4bcd,
  0x4bd5,0x4bdf,0x4bf1,0x4bfb,0x4c09,0x4c1d,0x4c3f,0x4c47,
  0x4c4b,0x4c55,0x4c7d,0x4c8d,0x4ca5,0x4cbd,0x4cc5,0x4cd1,
  0x4cd7,0x4ced,0x4d1f,0x4d31,0x4d45,0x4d4f,0x4d51,0x4d75,
  0x4d8f,0x4dad,0x4db3,0x4dc1,0x4dfb,0x4e01,0x4e15,0x4e37,
  0x4e43,0x4e49,0x4e51,0x4e6b,0x4e6d,0x4e83,0x4e9b,0x4ed3,
  0x4ee5,0x4ee9,0x4eef,0x4f1d,0x4f55,0x4f5f,0x4f69,0x4f6f,
  0x4f9f,0x4fa9,0x4fb7,0x4fbb,0x4feb,0x5007,0x501f,0x5025,
  0x505d,0x50a1,0x50d5,0x50e9,0x512d,0x5153,0x5159,0x515f,
  0x517d,0x518b,0x5199,0x51a5,0x51b1,0x51d1,0x51db,0x51eb,
  0x51ed,0x51f5,0x51ff,0x5205,0x520f,0x522d,0x5265,0x528b,
  0x5299,0x52a3,0x52b1,0x52c3,0x52c5,0x52d1,0x52e7,0x5325,
  0x5329,0x535b,0x539b,0x539d,0x53a7,0x53b9,0x53bf,0x53c7,
  0x53f1,0x5403,0x540f,0x5417,0x5435,0x5439,0x5447,0x545f,
  0x547b,0x548d,0x549f,0x54bb,0x54dd,0x54e1,0x5519,0x552f,
  0x5537,0x5579,0x557f,0x5583,0x5585,0x5591,0x55ad,0x55c1,
  0x55d3,0x55e9,0x55f7,0x560d,0x5615,0x5629,0x562f,0x5631,
  0x569b,0x56ab,0x56ad,0x56b3,0x56c1,0x56c7,0x56cd,0x56f7,
  0x572b,0x572d,0x573f,0x5759,0x5763,0x577d,0x5793,0x57af,
  0x57bd,0x57c3,0x57c5,0x57d7,0x5803,0x5821,0x5833,0x5835,
  0x583f,0x5841,0x5895,0x5899,0x58a3,0x58af,0x58dd,0x58e7,
  0x58eb,0x58ff,0x594f,0x597f,0x598f,0x5997,0x599d,0x59a1,
  0x59a7,0x59b5,0x59df,0x59e9,0x5a1f,0x5a25,0x5a3b,0x5a3d,
  0x5a43,0x5a45,0x5a49,0x5a75,0x5ab5,0x5ac1,0x5ad3,0x5ad5,
  0x5ad9,0x5b0f,0x5b17,0x5b2b,0x5b39,0x5b63,0x5b69,0x5b77,
  0x5b99,0x5bc3,0x5bc5,0x5be1,0x5beb,0x5c07,0x5c0b,0x5c19,
  0x5c31,0x5c49,0x5c5d,0x5c79,0x5c83,0x5ca1,0x5cbf,0x5cc1,
  0x5ccb,0x5ccd,0x5ce5,0x5d05,0x5d1d,0x5d33,0x5d6f,0x5d8d,
  0x5d95,0x5da3,0x5da9,0x5dbb,0x5dcf,0x5e03,0x5e0f,0x5e1b,
  0x5e2d,0x5e47,0x5e7d,0x5e81,0x5e99,0x5ed7,0x5edb,0x5ee7,
  0x5ef5,0x5f07,0x5f19,0x5f23,0x5f3d,0x5f43,0x5f45,0x5f61,
  0x5f6b,0x5f75,0x5f97,0x5fdf,0x600d,0x6015,0x602f,0x603d,
  0x6067,0x606b,0x6089,0x609d,0x60ab,0x60b3,0x60b9,0x60d5,
  0x6109,0x6111,0x612d,0x6139,0x6141,0x617d,0x619f,0x61a5,
  0x61af,0x61bb,0x61ed,0x61f5,0x6233,0x6263,0x626f,0x627d,
  0x628d,0x629f,0x62a5,0x62a9,0x62b7,0x62d7,0x62db,0x62dd,
  0x62f3,0x6323,0x6331,0x6337,0x635b,0x636d,0x6375,0x6389,
  0x6391,0x63a1,0x63bf,0x63ef,0x6409,0x6417,0x6441,0x6447,
  0x648b,0x649f,0x64e7,0x64f5,0x6501,0x650b,0x6537,0x653b,
  0x6545,0x655b,0x6591,0x65ad,0x65b9,0x65bf,0x65d5,0x65ef,
  0x65f7,0x6607,0x660d,0x6623,0x663b,0x665d,0x666b,0x6683,
  0x6697,0x66a7,0x66ab,0x66b5,0x66d9,0x66df,0x6711,0x671b,
  0x674b,0x675f,0x6769,0x6777,0x6781,0x67a3,0x67d1,0x6811,
  0x681d,0x6853,0x6877,0x6893,0x689f,0x68a5,0x68a9,0x68c5,
  0x68ff,0x6919,0x692f,0x696b,0x6973,0x698f,0x699d,0x69a1,
  0x69ab,0x69e9,0x6a01,0x6a2f,0x6a51,0x6a5b,0x6a67,0x6a6b,
  0x6a6d,0x6a75,0x6a83,0x6a97,0x6ab3,0x6ab5,0x6abf,0x6acb,
  0x6ae9,0x6aef,0x6b03,0x6b17,0x6b2b,0x6b33,0x6b39,0x6b47,
  0x6b4b,0x6b69,0x6b7d,0x6b81,0x6b8d,0x6bb7,0x6bc5,0x6bd7,
  0x6be1,0x6bed,0x6bf9,0x6c2f,0x6c3d,0x6c57,0x6c73,0x6cb5,
  0x6cb9,0x6cc1,0x6cc7,0x6ce5,0x6d09,0x6d17,0x6d2b,0x6d53,
  0x6d63,0x6d65,0x6d69,0x6da3,0x6db1,0x6dbb,0x6dbd,0x6dc5,
  0x6e09,0x6e27,0x6e2d,0x6e33,0x6e41,0x6e53,0x6e7b,0x6e81,
  0x6e87,0x6e95,0x6eaf,0x6ec3,0x6edb,0x6edd,0x6ef9,0x6f01,
  0x6f15,0x6f29,0x6f31,0x6f3b,0x6f4f,0x6f57,0x6f91,0x6fc7,
  0x6fd9,0x6fdf,0x6fe9,0x7005,0x701d,0x7033,0x707d,0x70a9,
  0x70bb,0x70c9,0x70cf,0x70e1,0x7113,0x7115,0x7119,0x7131,
  0x7137,0x713d,0x716b,0x719b,0x71a1,0x71b5,0x71df,0x71e5,
  0x71f7,0x71fb,0x721f,0x7231,0x723b,0x724f,0x7267,0x72b3,
  0x72c1,0x72cd,0x72df,0x72e5,0x72f1,0x72f7,0x7303,0x7309,
  0x7327,0x732b,0x738d,0x7393,0x73a5,0x73bd,0x73d1,0x73ff,
  0x7413,0x7415,0x744f,0x745b,0x746b,0x746d,0x74ab,0x74b3,
  0x74cd,0x74df,0x74e9,0x74fd,0x751b,0x7521,0x7577,0x757b,
  0x7599,0x75a3,0x75bd,0x75cf,0x75d1,0x75eb,0x75f5,0x7639,
  0x7647,0x7653,0x7655,0x765f,0x7663,0x7669,0x769f,0x76a3,
  0x76af,0x76d1,0x76eb,0x76f9,0x770b,0x7749,0x7757,0x776d,
  0x7773,0x77a7,0x77b5,0x77c7,0x77d3,0x77d5,0x7815,0x7825,
  0x7831,0x783d,0x786d,0x78cb,0x78cd,0x78d9,0x78ef,0x7917,
  0x7927,0x794d,0x7959,0x796f,0x7971,0x797b,0x7981,0x7987,
  0x79b1,0x79c9,0x79d7,0x79dd,0x7a03,0x7a1b,0x7a2b,0x7a35,
  0x7a3f,0x7a4b,0x7a55,0x7a81,0x7a8d,0x7ab7,0x7abb,0x7ac3,
  0x7ae1,0x7af5,0x7b23,0x7b4f,0x7b51,0x7b5d,0x7b79,0x7b8f,
  0x7ba1,0x7bab,0x7bb9,0x7bd3,0x7be3,0x7bf1,0x7c05,0x7c27,
  0x7c2d,0x7c59,0x7c8b,0x7c93,0x7c95,0x7ca3,0x7cb7,0x7cc3,
  0x7cd1,0x7cf9,0x7cff,0x7d01,0x7d15,0x7d37,0x7d45,0x7d73,
  0x7d79,0x7d91,0x7d97,0x7da7,0x7db3,0x7dc7,0x7dcd,0x7de9,
  0x7dfb,0x7dfd,0x7e0d,0x7e19,0x7e2f,0x7e61,0x7e75,0x7e7f,
  0x7e9d,0x7eab,0x7ed3,0x7ee3,0x7ee5,0x7f09,0x7f21,0x7f3f,
  0x7f4d,0x7f55,0x7f71,0x7f8b,0x7f8d,0x7f9f,0x7fc5,0x7fd1,
  0x7fe7,0x8003,0x8011,0x8017,0x802d,0x8035,0x805f,0x8077,
  0x8081,0x8087,0x8093,0x80a5,0x80c3,0x80cf,0x80dd,0x80e7,
  0x80f5,0x8101,0x8115,0x8125,0x8157,0x815d,0x8161,0x816d,
  0x8185,0x81a1,0x81a7,0x81b9,0x81cb,0x81cd,0x81df,0x81fd,
  0x8213,0x823b,0x8245,0x827f,0x8289,0x828f,0x829b,0x82cb,
  0x82d9,0x8317,0x831b,0x832b,0x8333,0x8347,0x834d,0x835f,
  0x8363,0x8369,0x8371,0x838b,0x8399,0x83af,0x83bd,0x83c5,
  0x83d1,0x8419,0x8423,0x842f,0x8431,0x8437,0x8467,0x846d,
  0x8479,0x8483,0x8497,0x84a1,0x84b5,0x84df,0x84f7,0x84fd,
  0x851d,0x8521,0x8527,0x8533,0x8547,0x854b,0x855f,0x8571,
  0x857b,0x8581,0x858d,0x85a3,0x85b1,0x85c5,0x85c9,0x85db,
  0x85ed,0x85f3,0x8609,0x8611,0x861d,0x862b,0x8655,0x8659,
  0x8665,0x867d,0x8681,0x86a9,0x86af,0x86b7,0x86c3,0x86db,
  0x86ff,0x870b,0x870d,0x8715,0x8729,0x8731,0x876b,0x8785,
  0x8789,0x878f,0x87ab,0x87ad,0x87b3,0x87bf,0x87df,0x87f1,
  0x8801,0x8837,0x8849,0x885d,0x8861,0x8879,0x8883,0x888f,
  0x88b5,0x88c7,0x88d9,0x8905,0x8909,0x8917,0x8927,0x8965,
  0x8999,0x89af,0x89cf,0x89d1,0x89dd,0x89ed,0x89ff,0x8a39,
  0x8a55,0x8a63,0x8a65,0x8a71,0x8a7b,0x8a87,0x8a8b,0x8a9f,
  0x8aaf,0x8ab1,0x8add,0x8b0b,0x8b0d,0x8b3b,0x8b43,0x8b6d,
  0x8b73,0x8b91,0x8b9d,0x8bc1,0x8bcb,0x8bd5,0x8be5,0x8bf1,
  0x8c21,0x8c27,0x8c47,0x8c55,0x8c63,0x8c69,0x8c87,0x8c93,
  0x8ca3,0x8ca9,0x8cbb,0x8cc9,0x8ce1,0x8d07,0x8d0b,0x8d25,
  0x8d43,0x8d51,0x8d5b,0x8d7f,0x8d85,0x8d89,0x8d9b,0x8da1,
  0x8dbf,0x8dcb,0x8dcd,0x8de5,0x8de9,0x8def,0x8e0b,0x8e1f,
  0x8e3b,0x8e51,0x8e57,0x8e73,0x8e75,0x8e8f,0x8e97,0x8e9b,
  0x8ea1,0x8ebf,0x8ec1,0x8ed3,0x8ee5,0x8efb,0x8f0f,0x8f27,
  0x8f2b,0x8f39,0x8f65,0x8f77,0x8f7d,0x8fa5,0x8fc9,0x8fd1,
  0x8fe1,0x900b,0x9019,0x901f,0x9025,0x9037,0x903d,0x9043,
  0x9057,0x9061,0x906d,0x9075,0x9089,0x9091,0x909b,0x90ad,
  0x9105,0x9109,0x9117,0x911b,0x911d,0x9141,0x914d,0x917b,
  0x9199,0x91a3,0x91a9,0x91b1,0x91bd,0x91cf,0x91e1,0x9203,
  0x9211,0x921b,0x9227,0x9233,0x923f,0x924d,0x9253,0x927b,
  0x9287,0x9295,0x92a3,0x92b7,0x92c9,0x92dd,0x92f5,0x92f9,
  0x9313,0x9331,0x9345,0x9349,0x9357,0x9373,0x937f,0x9383,
  0x939d,0x93a1,0x93b5,0x93cb,0x93df,0x93ef,0x93f1,0x93fb,
  0x940f,0x942b,0x9435,0x944b,0x9453,0x947d,0x9493,0x9499,
  0x94a5,0x94a9,0x94bb,0x94d7,0x94db,0x94e1,0x94e7,0x94f5,
  0x94ff,0x9507,0x9515,0x9529,0x9531,0x953d,0x9545,0x9561,
  0x956b,0x9589,0x959b,0x95ad,0x95df,0x95e9,0x9613,0x9615,
  0x9631,0x9637,0x963d,0x9643,0x9667,0x9679,0x9697,0x96ad,
  0x96c1,0x96cb,0x96df,0x96fb,0x972b,0x9735,0x974d,0x9763,
  0x9795,0x979f,0x97a3,0x97a5,0x97a9,0x97b1,0x97db,0x97dd,
  0x97f3,0x97ff,0x9809,0x9821,0x9847,0x9895,0x9899,0x98af,
  0x98bb,0x98c5,0x98d7,0x98dd,0x990b,0x9919,0x9925,0x9929,
  0x992f,0x9937,0x9943,0x9945,0x996d,0x997f,0x9989,0x9991,
  0x999b,0x99b9,0x99bf,0x99c1,0x99cb,0x99d5,0x99e3,0x9a25,
  0x9a4f,0x9a57,0x9a5b,0x9a5d,0x9a61,0x9a6d,0x9a79,0x9a7f,
  0x9a97,0x9a9d,0x9aab,0x9ab9,0x9ac7,0x9ae3,0x9afb,0x9b11,
  0x9b1b,0x9b33,0x9b41,0x9b53,0x9b5f,0x9b6f,0x9b77,0x9b8b,
  0x9b93,0x9bb7,0x9bbd,0x9beb,0x9bff,0x9c15,0x9c2f,0x9c3b,
  0x9c51,0x9c5d,0x9c67,0x9c75,0x9c79,0x9c85,0x9cbf,0x9cd5,
  0x9ce5,0x9cf1,0x9cf7,0x9d03,0x9d55,0x9d59,0x9d63,0x9d81,
  0x9d8b,0x9d99,0x9d9f,0x9da5,0x9dbd,0x9e03,0x9e11,0x9e17,
  0x9e21,0x9e39,0x9e47,0x9e55,0x9e59,0x9e63,0x9e69,0x9e9f,
  0x9ea3,0x9eed,0x9ef3,0x9f0b,0x9f13,0x9f23,0x9f2f,0x9f49,
  0x9f4f,0x9f67,0x9f9b,0x9fad,0x9fb3,0x9fb5,0x9fcd,0xa00b,
  0xa023,0xa025,0xa02f,0xa043,0xa05b,0xa05d,0xa067,0xa06b,
  0xa06d,0xa089,0xa091,0xa0b5,0xa0cd,0xa0df,0xa0fb,0xa103,
  0xa10f,0xa11d,0xa127,0xa133,0xa135,0xa139,0xa147,0xa16f,
  0xa17b,0xa181,0xa19f,0xa1b1,0xa1c5,0xa1cf,0xa1e1,0xa20f,
  0xa241,0xa247,0xa263,0xa26f,0xa299,0xa2a9,0xa2c5,0xa2c9,
  0xa2d7,0xa2db,0xa2ed,0xa307,0xa30b,0xa30d,0xa319,0xa323,
  0xa345,0xa36b,0xa373,0xa385,0xa38f,0xa397,0xa3a1,0xa3b3,
  0xa3c1,0xa3cd,0xa3d5,0xa403,0xa405,0xa409,0xa453,0xa459,
  0xa481,0xa48b,0xa499,0xa4a5,0xa4af,0xa4b1,0xa4bd,0xa4c3,
  0xa4db,0xa4ed,0xa501,0xa523,0xa525,0xa529,0xa54f,0xa583,
  0xa5a7,0xa5b9,0xa5cd,0xa5e9,0xa5f1,0xa613,0xa62f,0xa637,
  0xa651,0xa661,0xa673,0xa691,0xa69b,0xa69d,0xa6c7,0xa6f1,
  0xa6fd,0xa71b,0xa72b,0xa739,0xa75f,0xa771,0xa78d,0xa7b1,
  0xa7bb,0xa7bd,0xa7cf,0xa7d1,0xa7dd,0xa7e7,0xa82b,0xa833,
  0xa839,0xa84b,0xa84d,0xa855,0xa869,0xa87d,0xa881,0xa893,
  0xa89f,0xa8a3,0xa8a9,0xa8b7,0xa8cf,0xa8dd,0xa8e1,0xa90b,
  0xa913,0xa919,0xa949,0xa975,0xa9a7,0xa9ad,0xa9e9,0xa9ef,
  0xa9f7,0xaa0d,0xaa15,0xaa2f,0xaa31,0xaa43,0xaa4f,0xaa51,
  0xaa61,0xaa79,0xaa9d,0xaaad,0xaab9,0xaad3,0xab03,0xab0f,
  0xab1d,0xab35,0xab39,0xab47,0xab63,0xab8b,0xab99,0xabbb,
  0xabc5,0xabd1,0xabd7,0xabdb,0xabf3,0xac01,0xac07,0xac0d,
  0xac29,0xac37,0xac5b,0xac6d,0xac73,0xac85,0xac8f,0xac97,
  0xacc7,0xacd5,0xacdf,0xace9,0xacfd,0xad05,0xad11,0xad1b,
  0xad21,0xad6f,0xada3,0xadaf,0xadc9,0xaddd,0xade7,0xadeb,
  0xadf9,0xae09,0xae1d,0xae2d,0xae39,0xae3f,0xae5f,0xae63,
  0xae71,0xae77,0xae8d,0xae95,0xaea3,0xaebb,0xaed7,0xaedd,
  0xaee7,0xaef5,0xaf01,0xaf13,0xaf1f,0xaf29,0xaf2f,0xaf3b,
  0xaf43,0xaf49,0xaf67,0xaf75,0xaf8f,0xafab,0xafd3,0xafe3,
  0xaffd,0xb00f,0xb01b,0xb02d,0xb035,0xb03f,0xb055,0xb077,
  0xb093,0xb0af,0xb0c3,0xb0c5,0xb0d1,0xb0e1,0xb0f3,0xb107,
  0xb10b,0xb13d,0xb157,0xb15d,0xb175,0xb18f,0xb197,0xb19d,
  0xb1a1,0xb1c7,0xb1df,0xb1e5,0xb20b,0xb213,0xb215,0xb243,
  0xb249,0xb24f,0xb27f,0xb289,0xb29b,0xb2c1,0xb2df,0xb2e9,
  0xb2ef,0xb2f7,0xb305,0xb317,0xb327,0xb37b,0xb381,0xb393,
  0xb3a5,0xb3b1,0xb3b7,0xb3c5,0xb3d7,0xb3db,0xb3dd,0xb3f9,
  0xb401,0xb40d,0xb443,0xb457,0xb46d,0xb475,0xb47f,0xb4a7,
  0xb4bf,0xb4ef,0xb4f7,0xb509,0xb51d,0xb53f,0xb547,0xb553,
  0xb555,0xb569,0xb577,0xb57b,0xb595,0xb5a9,0xb5bb,0xb5c3,
  0xb5d7,0xb5e1,0xb5e7,0xb5eb,0xb5f9,0xb605,0xb609,0xb617,
  0xb61b,0xb621,0xb627,0xb62d,0xb633,0xb635,0xb64b,0xb659,
  0xb67b,0xb681,0xb68b,0xb699,0xb6af,0xb6d1,0xb6e7,0xb6eb,
  0xb71f,0xb725,0xb745,0xb757,0xb76b,0xb779,0xb783,0xb791,
  0xb79d,0xb7a1,0xb7b3,0xb7fd,0xb84f,0xb861,0xb875,0xb885,
  0xb889,0xb89d,0xb8a1,0xb8ad,0xb8bf,0xb8d3,0xb8d5,0xb8fb,
  0xb903,0xb91d,0xb927,0xb92b,0xb947,0xb955,0xb959,0xb965,
  0xb977,0xb98d,0xb9af,0xb9bb,0xb9c3,0xb9c9,0xb9d1,0xb9ed,
  0xb9f3,0xba05,0xba11,0xba39,0xba59,0xba81,0xba87,0xba8d,
  0xba93,0xbab7,0xbad7,0xbb01,0xbb15,0xbb19,0xbb23,0xbb3d,
  0xbb43,0xbb49,0xbb51,0xbb5b,0xbb67,0xbb75,0xbb91,0xbbb5,
  0xbbbf,0xbbcb,0xbbcd,0xbbe5,0xbbe9,0xbc09,0xbc33,0xbc69,
  0xbc7b,0xbc8b,0xbc8d,0xbc9f,0xbca9,0xbcc3,0xbcdd,0xbd13,
  0xbd25,0xbd2f,0xbd67,0xbd7f,0xbd83,0xbd89,0xbdb9,0xbdc1,
  0xbdd3,0xbdd9,0xbde5,0xbe15,0xbe23,0xbe29,0xbe37,0xbe57,
  0xbe5b,0xbe61,0xbe6b,0xbe7f,0xbe83,0xbeab,0xbef1,0xbf21,
  0xbf27,0xbf33,0xbf35,0xbf65,0xbf81,0xbf8b,0xbfa3,0xbfb7,
  0xbfc3,0xbfed,0xbff5,0xbfff,0xc001,0xc007,0xc013,0xc01f,
  0xc025,0xc049,0xc073,0xc079,0xc085,0xc097,0xc09d,0xc0a7,
  0xc0b9,0xc0cb,0xc0d5,0xc0e3,0xc0fb,0xc111,0xc121,0xc133,
  0xc13f,0xc15f,0xc177,0xc17b,0xc17d,0xc193,0xc1a5,0xc1bd,
  0xc1c3,0xc1c9,0xc1db,0xc1e7,0xc1ed,0xc205,0xc209,0xc217,
  0xc22d,0xc24d,0xc255,0xc269,0xc299,0xc2b1,0xc2d1,0xc2dd,
  0xc2f5,0xc301,0xc30d,0xc325,0xc33d,0xc357,0xc35b,0xc361,
  0xc373,0xc383,0xc39b,0xc39d,0xc3ad,0xc3c7,0xc3e3,0xc3ef,
  0xc3fd,0xc405,0xc421,0xc433,0xc447,0xc44b,0xc45f,0xc46f,
  0xc47d,0xc493,0xc49f,0xc4a5,0xc4b7,0xc4c5,0xc4dd,0xc4eb,
  0xc4f9,0xc4ff,0xc50b,0xc515,0xc52f,0xc531,0xc549,0xc55b,
  0xc575,0xc579,0xc589,0xc59b,0xc5a1,0xc5b3,0xc5b5,0xc5bf,
  0xc5df,0xc5e9,0xc5ef,0xc5fd,0xc631,0xc645,0xc651,0xc66b,
  0xc675,0xc679,0xc68f,0xc697,0xc6ab,0xc6b9,0xc6c1,0xc6c7,
  0xc6d5,0xc6e9,0xc703,0xc753,0xc759,0xc777,0xc799,0xc79f,
  0xc7c3,0xc7d7,0xc7db,0xc7f5,0xc803,0xc81b,0xc841,0xc84d,
  0xc85f,0xc865,0xc869,0xc895,0xc8bd,0xc8c9,0xc8f5,0xc8f9,
  0xc8ff,0xc901,0xc90d,0xc915,0xc91f,0xc923,0xc929,0xc931,
  0xc937,0xc93b,0xc95b,0xc95d,0xc983,0xc9c7,0xc9cb,0xc9cd,
  0xc9d9,0xca25,0xca29,0xca2f,0xca49,0xca67,0xca8f,0xcaad,
  0xcabf,0xcad9,0xcae3,0xcaef,0xcb17,0xcb1b,0xcb1d,0xcb2b,
  0xcb3f,0xcb55,0xcb71,0xcb7b,0xcb87,0xcbbd,0xcbe7,0xcbf5,
  0xcc15,0xcc23,0xcc3d,0xcc49,0xcc6d,0xcc7f,0xcc83,0xcc85,
  0xcc9b,0xcca1,0xcca7,0xccb3,0xccc1,0xccd9,0xccf7,0xccfb,
  0xccfd,0xcd2b,0xcd33,0xcd87,0xcd9f,0xcda3,0xcdc5,0xcdcf,
  0xcdd7,0xcde1,0xcdeb,0xcded,0xcdf9,0xcdff,0xce03,0xce27,
  0xce35,0xce65,0xce71,0xce9f,0xcebb,0xcec3,0xcec5,0xcec9,
  0xced1,0xced7,0xcef3,0xceff,0xcf0d,0xcf1f,0xcf2f,0xcf73,
  0xcf79,0xcf8f,0xcf9d,0xcfa1,0xcfc7,0xcfd5,0xcfe9,0xcff7,
  0xd005,0xd009,0xd03f,0xd04d,0xd071,0xd077,0xd087,0xd08b,
  0xd08d,0xd095,0xd099,0xd0a3,0xd0b1,0xd0c5,0xd0cf,0xd0d1,
  0xd0d7,0xd0e1,0xd0f9,0xd10b,0xd125,0xd12f,0xd13d,0xd151,
  0xd15b,0xd16d,0xd19b,0xd1b9,0xd1bf,0xd1c1,0xd1d5,0xd1d9,
  0xd1ef,0xd1fd,0xd207,0xd215,0xd223,0xd229,0xd237,0xd257,
  0xd26b,0xd26d,0xd28f,0xd2a1,0xd303,0xd317,0xd327,0xd341,
  0xd35f,0xd369,0xd377,0xd37b,0xd381,0xd393,0xd399,0xd3b1,
  0xd3c9,0xd3d1,0xd3dd,0xd3e7,0xd415,0xd429,0xd43b,0xd461,
  0xd46b,0xd497,0xd49d,0xd4ab,0xd4b3,0xd4bf,0xd4c1,0xd4d3,
  0xd4e5,0xd4e9,0xd4f1,0xd50f,0xd527,0xd52b,0xd547,0xd559,
  0xd563,0xd57d,0xd5b7,0xd5e1,0xd5e7,0xd5f5,0xd605,0xd627,
  0xd62b,0xd64b,0xd663,0xd67d,0xd6a9,0xd6b7,0xd6c5,0xd6d7,
  0xd6e1,0xd6ed,0xd723,0xd76d,0xd77f,0xd7ad,0xd7b3,0xd7b5,
  0xd7bf,0xd7d9,0xd7f7,0xd7fb,0xd80d,0xd813,0xd849,0xd857,
  0xd86d,0xd889,0xd88f,0xd89b,0xd8a7,0xd8b5,0xd8c1,0xd8d3,
  0xd8d9,0xd8e5,0xd909,0xd91b,0xd933,0xd941,0xd947,0xd94d,
  0xd95f,0xd965,0xd971,0xd98b,0xd999,0xd99f,0xd9a3,0xd9a9,
  0xd9b1,0xd9c3,0xd9d7,0xd9f9,0xda05,0xda17,0xda27,0xda35,
  0xda3f,0xda59,0xda7d,0xda8b,0xda93,0xdaa3,0xdab1,0xdac3,
  0xdadd,0xdb1f,0xdb25,0xdb29,0xdb3b,0xdb45,0xdb4f,0xdb61,
  0xdb83,0xdba1,0xdbc7,0xdbcd,0xdbd5,0xdbdf,0xdbe3,0xdbe9,
  0xdbef,0xdc2b,0xdc39,0xdc41,0xdc6f,0xdc71,0xdc93,0xdc9f,
  0xdcaf,0xdcd1,0xdcd7,0xdcdb,0xdcf5,0xdd19,0xdd29,0xdd31,
  0xdd57,0xdd67,0xdd73,0xdd75,0xdd9d,0xddad,0xddd5,0xdde5,
  0xddf7,0xddfb,0xde07,0xde2f,0xde3d,0xde49,0xde4f,0xde51,
  0xde6d,0xde83,0xde85,0xde89,0xdea1,0xdead,0xdecb,0xdecd,
  0xded3,0xdf03,0xdf05,0xdf17,0xdf1d,0xdf33,0xdf59,0xdf69,
  0xdf6f,0xdf71,0xdfbb,0xdfc9,0xdfeb,0xe003,0xe00f,0xe035,
  0xe047,0xe04b,0xe05f,0xe07b,0xe08d,0xe0a9,0xe0b1,0xe0b7,
  0xe0c5,0xe0d7,0xe101,0xe10b,0xe131,0xe149,0xe14f,0xe151,
  0xe15d,0xe167,0xe17f,0xe1b3,0xe1d3,0xe1ef,0xe207,0xe219,
  0xe223,0xe231,0xe245,0xe24f,0xe267,0xe279,0xe285,0xe29b,
  0xe29d,0xe2a1,0xe2ab,0xe2ad,0xe2bf,0xe2c1,0xe2d5,0xe30f,
  0xe311,0xe335,0xe359,0xe363,0xe365,0xe377,0xe38d,0xe393,
  0xe39f,0xe3af,0xe3b7,0xe3c3,0xe3db,0xe3f3,0xe3ff,0xe41f,
  0xe431,0xe449,0xe45b,0xe46b,0xe46d,0xe473,0xe485,0xe491,
  0xe497,0xe49d,0xe4a1,0xe4ab,0xe4cb,0xe4cd,0xe4f1,0xe4fd,
  0xe503,0xe517,0xe51b,0xe52d,0xe533,0xe581,0xe595,0xe5a5,
  0xe605,0xe621,0xe639,0xe63f,0xe647,0xe653,0xe669,0xe687,
  0xe6bb,0xe6bd,0xe6d7,0xe6dd,0xe6f5,0xe6f9,0xe701,0xe729,
  0xe76d,0xe775,0xe783,0xe78f,0xe7ab,0xe7ad,0xe7b5,0xe7cb,
  0xe7d3,0xe7e5,0xe7f7,0xe801,0xe82f,0xe843,0xe85b,0xe86d,
  0xe879,0xe889,0xe891,0xe8a7,0xe8bf,0xe8c1,0xe8cb,0xe8cd,
  0xe8d3,0xe8fb,0xe903,0xe90f,0xe921,0xe927,0xe92b,0xe935,
  0xe959,0xe95f,0xe963,0xe969,0xe971,0xe98d,0xe9b7,0xe9c5,
  0xe9cf,0xe9d7,0xea09,0xea1b,0xea2d,0xea4b,0xea59,0xea71,
  0xea7d,0xea81,0xea8d,0xeabb,0xeac3,0xeac9,0xeaed,0xeb07,
  0xeb0b,0xeb19,0xeb29,0xeb3b,0xeb45,0xeb5d,0xeb67,0xeb6b,
  0xeb73,0xeb75,0xeb97,0xeb9b,0xebad,0xebb3,0xebcd,0xebd5,
  0xebe3,0xec09,0xec11,0xec21,0xec35,0xec4b,0xec65,0xec69,
  0xec7d,0xec93,0xec99,0xecb7,0xeccf,0xecff,0xed07,0xed15,
  0xed23,0xed37,0xed49,0xed5d,0xed61,0xed6b,0xed97,0xedab,
  0xedc7,0xedcd,0xedd9,0xeddf,0xedef,0xedfd,0xee01,0xee0b,
  0xee0d,0xee1f,0xee75,0xee83,0xee9d,0xeead,0xeec7,0xeecb,
  0xeed9,0xeee3,0xeef1,0xeef7,0xef0f,0xef21,0xef2d,0xef33,
  0xef39,0xef4d,0xef77,0xef95,0xef9f,0xefbb,0xefe7,0xefeb,
  0xeff3,0xf007,0xf00d,0xf01f,0xf029,0xf02f,0xf045,0xf085,
  0xf097,0xf0ab,0xf0bf,0xf0c7,0xf0d5,0xf0df,0xf0f1,0xf0f7,
  0xf111,0xf11b,0xf135,0xf141,0xf14b,0xf153,0xf163,0xf171,
  0xf18d,0xf1c5,0xf1e1,0xf1e7,0xf1f3,0xf1f5,0xf21d,0xf247,
  0xf24d,0xf255,0xf259,0xf26f,0xf27b,0xf287,0xf29f,0xf2a5,
  0xf2db,0xf2f9,0xf2ff,0xf301,0xf30b,0xf315,0xf32f,0xf337,
  0xf385,0xf389,0xf391,0xf397,0xf3b3,0xf3df,0xf3e5,0xf405,
  0xf40f,0xf417,0xf421,0xf439,0xf453,0xf455,0xf465,0xf47b,
  0xf48b,0xf499,0xf4a3,0xf4bd,0xf4cf,0xf4f3,0xf4f5,0xf4f9,
  0xf50d,0xf519,0xf525,0xf53b,0xf551,0xf561,0xf56d,0xf591,
  0xf59d,0xf5b5,0xf5bf,0xf5c1,0xf5c7,0xf61f,0xf623,0xf63b,
  0xf645,0xf64f,0xf685,0xf6b5,0xf6d9,0xf6ef,0xf6fb,0xf72d,
  0xf73f,0xf74d,0xf753,0xf76f,0xf787,0xf78b,0xf795,0xf7a3,
  0xf7b1,0xf7b7,0xf7c3,0xf7c9,0xf7db,0xf7ff,0xf803,0xf809,
  0xf80f,0xf827,0xf83f,0xf86f,0xf871,0xf877,0xf893,0xf8db,
  0xf8ed,0xf8f3,0xf8f5,0xf915,0xf923,0xf93b,0xf93d,0xf94f,
  0xf951,0xf973,0xf979,0xf985,0xf99b,0xf9b3,0xf9b9,0xf9c7,
  0xf9e3,0xf9e9,0xf9f7,0xfa01,0xfa07,0xfa13,0xfa23,0xfa75,
  0xfa83,0xfa97,0xfa9b,0xfaa1,0xfac1,0xfacb,0xfad9,0xfadf,
  0xfae5,0xfb05,0xfb0f,0xfb21,0xfb35,0xfb4d,0xfb5f,0xfb69,
  0xfb81,0xfb8d,0xfba3,0xfba9,0xfbb7,0xfbc9,0xfbcf,0xfbdb,
  0xfbe1,0xfc0b,0xfc0d,0xfc1f,0xfc49,0xfc5b,0xfc67,0xfc75,
  0xfc83,0xfcad,0xfcd3,0xfcef,0xfd0f,0xfd17,0xfd1d,0xfd2b,
  0xfd2d,0xfd39,0xfd47,0xfd53,0xfd71,0xfd8b,0xfd99,0xfda3,
  0xfdaf,0xfdb1,0xfddd,0xfde1,0xfdeb,0xfe2d,0xfe33,0xfe41,
  0xfe4d,0xfe59,0xfe7d,0xfe87,0xfe99,0xfeb1,0xfebd,0xfec9,
  0xfeeb,0xfeff,0xff13,0xff23,0xff29,0xff37,0xff4f,0xff61,
  0xff73,0xff7f,0xff91,0xffb3,0xffc7,0xffd9,0xffe9,0xffef,
  0xfffd};
