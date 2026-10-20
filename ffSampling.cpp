
static int convert_case(int lvl)
{
	return lvl << (lvl - 1);
}

static int convert(int j, int lvl, int z)
{
	if (lvl == 0) {
		return 1;
	}
	z ^= 1;
	int out = 0;
	// Keep MSB (bit lvl)
	out |= j & (1 << lvl);
	// Set bit lvl-1 to z
	out |= (z & 1) << (lvl - 1);
	// Fill the remaining bits:
	// out[i] = ~j[i+1], for i = 0 ... lvl-2
	for (int i = 0; i <= lvl - 2; i++) {
		int bit = (((j >> (i + 1)) & 1) ^ 1);
		out |= bit << i;
	}
	return out;
}
static int convert2(int j, int lvl)
{
	int p, y;
	p = 9 - lvl; //	printf("p  %d  ", p);
	y = 0;
	int t = 0;
	for (int i = p - 1; i >= 0; i--) {
		if ((j >> i) & 1) {
			y += convert_case(9 - t);
		}
		t++;
	}
//	printf("\n");
	return y;
}

static int convert3(int lvl)
{
	int y = 0, p = 9 - lvl;
	if (p == 0) y = 0;
	else if (p == 1) y = 512;   // 256
	else if (p == 2) y = 768;   // 128
	else if (p == 3) y = 896;   //  64
	else if (p == 4) y = 960;   //  32
	else if (p == 5) y = 992;   //  16
	else if (p == 6) y = 1008;  //  8
	else if (p == 7) y = 1016;
	else if (p == 8) y = 1020;
	return y;
}


static int convert_fft(int j, int lvl)
{
	int out = 0;
	for (int i = 0; i < lvl; i++) {
		int bit = ((j >> i) & 1) ^ 1;     // Invert bit i
		out |= bit << (lvl - 1 - i);      // Place in reversed position
	}
	return out;
}

void ffSamplerZ(samplerZ samp, void* samp_ctx, fpr tx[1024], fpr z[512], const fpr tree, int i) {
	fpr t[2];
	fpr sigma;
	int index = i;
	t[0] = z[index]; t[1] = z[index + 256]; 
	printf("i %d  index %d  tree -> 3 %lf  ", i, index, tree);
	printf("t1 -> t2 %lf  %lf  \n", t[0], t[1]);
	sigma = tree;
	SamplePair r;
	r = sampler2(samp_ctx, t[0], t[1], sigma);
	z[index] = fpr_of(r.s0);
	z[index+256] = fpr_of(r.s1);
	printf("z1 -> z3 %lf  %lf \n", z[index], z[index + 256]);
}

static void ffSampling(samplerZ samp, void* samp_ctx, fpr z0[512], fpr z1[512], const fpr* tree, const fpr* t0, const fpr* t1) {
	size_t n, lvl, index, index1, idx = 0, p = 0;
	int state = 0;
	int i1 = 0, i2 = 0, i3 = 0;
	fpr tree0[512], tree1[512], treex[512], tree2;
	fpr tx[512];
	n = 512;
	lvl = 9;
	for (int i = 0; i < n; i++) {
		z1[i] = t1[i];
		z0[i] = t1[i];
		treex[i] = fpr_of(0);
	}
	//	for (int i = 0; i < 256; i++) {
	//		printf("z1[%d]  %lf  ", i, z1[i]);
	//	}printf("\n\n");
	idx = 1;
	while (1) {
		if (state == 0) {
			n = 1 << lvl;
			index = convert(idx, 9 - lvl,0); index1 = convert_fft(idx, 9 - lvl);
			poly_split_fft9(tx,z1, index1, lvl);
			lvl = lvl - 1;
			idx = idx << 1;
			if (lvl == 1) state = 1;
		}
		if (state == 1) {
			n = 1 << lvl;
			index = convert(idx, 9 - lvl,0);
			index1 = convert_fft(idx, 9 - lvl);
			i1 = convert3(lvl) + convert2(index, lvl); i2 = i1 + n; i3 = i1 + n + (lvl << (lvl - 1));
			if (idx % 2 == 0) tree2 = tree[i3];
			else tree2 = tree[i2];
			ffSamplerZ(samp, samp_ctx,tx, z1, tree2, index1);
			if (idx % 2 == 0) {
				state = 2;
			}
			else {
				idx = idx >> 1; lvl++;
				state = 3;
			}
		}
		if (state == 2) {
			n = 1 << lvl; p = 1 << (9 - lvl);
			index1 = convert_fft(idx, 9 - lvl);
			poly_sub_logn(tx, z1, lvl, index1);
			idx++;
			index = convert(idx, 9 - lvl, 0);
			i1 = convert3(lvl) + convert2(index, lvl); i2 = i1 + n;
			index1 = convert_fft(idx, 9 - lvl);
			for (int i = 0; i < n; i++) {
				treex[i * p + index1] = tree[i1 + i];
			}
			poly_mult_logn(tx, treex, lvl, index1);
			poly_add_logn(tx, z1, lvl, index1);
			if (lvl != 1) state = 0;
			else {
			  state = 1;
			}
		}
		if (state == 3) {
			index = convert(idx, 9 - lvl,0); index1 = convert_fft(idx, 9 - lvl);
			poly_merge_fft9(z1, index1, lvl);
			if (lvl == 9) break;
			if (idx % 2 == 1) {
				state = 3;
				lvl++; idx = idx >> 1;
			}
			else {
				state = 2; 
			}
		}
	}
	n = 1 << lvl; p = 1 << (9 - lvl);
	index = convert(idx, 9 - lvl,1);  i1 = convert3(lvl) + convert2(index, lvl); i2 = i1 + n; i3 = i1 + n + (lvl << (lvl - 1));
	index1 = convert_fft(idx, 9 - lvl);
	poly_sub(z0, z1, lvl);
	for (int i = 0; i < n; i++) {
		treex[i * p + index1] = tree[i1 + i];
	}
	poly_mult_fft(z0, treex, lvl);
	poly_add(z0, t0, lvl);
	idx = 1;state = 0;
	while (1) {
		if (state == 0) {
			n = 1 << lvl;
			index = convert(idx, 9 - lvl, 1); index1 = convert_fft(idx, 9 - lvl);
			i1 = convert3(lvl) + convert2(index, lvl); 
			poly_split_fft9(tx, z0, index1, lvl);
			lvl = lvl - 1;
			idx = idx << 1;
			if (lvl == 1) state = 1;
		}
		if (state == 1) {
			n = 1 << lvl;
			index = convert(idx, 9 - lvl, 1);
			index1 = convert_fft(idx, 9 - lvl);
			i1 = convert3(lvl) + convert2(index, lvl); i2 = i1 + n; i3 = i1 + n + (lvl << (lvl - 1));
			if (idx % 2 == 0) tree2 = tree[i3];
			else tree2 = tree[i2];
			ffSamplerZ(samp, samp_ctx, tx, z0, tree2, index1);
			if (idx % 2 == 0) {
				state = 2;
			}
			else {
				idx = idx >> 1; lvl++;
				state = 3;
			}
		}
		if (state == 2) {
			n = 1 << lvl; p = 1 << (9 - lvl);
			index1 = convert_fft(idx, 9 - lvl);
			poly_sub_logn(tx, z0, lvl, index1);
			idx++;
			index = convert(idx, 9 - lvl, 1);
			i1 = convert3(lvl) + convert2(index, lvl); i2 = i1 + n;
			index1 = convert_fft(idx, 9 - lvl);
			for (int i = 0; i < n; i++) {
				treex[i * p + index1] = tree[i1 + i];
			}
			poly_mult_logn(tx, treex, lvl, index1);
			poly_add_logn(tx, z0, lvl, index1);
			if (lvl != 1) state = 0;
			else {
				state = 1;
			}
		}
		if (state == 3) {
			index = convert(idx, 9 - lvl, 1); index1 = convert_fft(idx, 9 - lvl);
			poly_merge_fft9(z0, index1, lvl);
			if (lvl == 9) break;
			if (idx % 2 == 1) {
				state = 3;
				lvl++; idx = idx >> 1;
			}
			else {
				state = 2;
			}
		}
	}
	return;
}
