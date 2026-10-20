
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
		 //	i1 = convert3(lvl) + convert2(index, lvl); // i3 = i1 + n + (lvl << (lvl - 1));
			printf("state 0 convert 3 %d  convert2 %d   state %d  ", convert3(lvl), convert2(index, lvl), state);
			printf("j %d index1 %d  index  %d   lvl  %d   n % d   p %d   \n", idx, index1, index, lvl, n, p);
			//poly_copy_fft9(tx, z1, lvl, index);
			poly_split_fft9(tx,z1, index1, lvl);
			//poly_copy_fft9(ty, z1, lvl, j);
			lvl = lvl - 1;
			idx = idx << 1;
			if (lvl == 1) state = 1;
		}
		if (state == 1) {
			n = 1 << lvl;
			index = convert(idx, 9 - lvl,0);
			index1 = convert_fft(idx, 9 - lvl);
			i1 = convert3(lvl) + convert2(index, lvl); i2 = i1 + n; i3 = i1 + n + (lvl << (lvl - 1));
			printf("\n state1 index %d  j %d  index1 %d   lvl %d  i1 %d  i2 %d  i3 %d   %d   %d\n", index, idx,index1, lvl,i1,i2,i3, convert3(lvl), convert2(index, lvl));
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
			//printf(" \n");
			poly_sub_logn(tx, z1, lvl, index1);
			idx++;
			index = convert(idx, 9 - lvl, 0);
			i1 = convert3(lvl) + convert2(index, lvl); i2 = i1 + n;
			index1 = convert_fft(idx, 9 - lvl);
			printf("\n state2  p %d  index %d  j %d index1 %d  lvl %d  i1 %d  i2 %d  conv3 %d \n", p, index, idx, index1, lvl, i1, i2, convert3(lvl));
			for (int i = 0; i < n; i++) {
				treex[i * p + index1] = tree[i1 + i];
				//printf(" treex[%d]  %lf  ,", i * p + index1,treex[i*p + index1]);
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
			//printf("\n"); 
			printf("\n state3  index %d  j %d  index1 %d  lvl %d  i1 %d  i2 %d  i3 %d  \n", index, idx, index1, lvl, i1, i2, i3);
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
	printf("j %d index  %d  index1  %d lvl  %d   n % d   p %d   i1  %d   i2 %d i3 %d\n", idx,  index, index1, lvl, n, p, i1, i2, i3);
	poly_sub(z0, z1, lvl);
	for (int i = 0; i < n; i++) {
		treex[i * p + index1] = tree[i1 + i];
	}
	poly_mult_fft(z0, treex, lvl);
	poly_add(z0, t0, lvl);
	for (int i = 0; i < 512; i++) {
		//printf("z0[%d]  %lf  ", i, z0[i]);
	}printf("\n\n");
	idx = 1;state = 0;
	while (1) {
		if (state == 0) {
			n = 1 << lvl;
			index = convert(idx, 9 - lvl, 1); index1 = convert_fft(idx, 9 - lvl);
			i1 = convert3(lvl) + convert2(index, lvl); 
			printf("state 0 convert 3 %d  convert2 %d   state %d  ", convert3(lvl), convert2(index, lvl), state);
			printf("j %d index1 %d  index  %d   lvl  %d   n % d   p %d   i1  %d   \n", idx, index1, index, lvl, n, p, i1);
			//poly_copy_fft9(tx, z1, lvl, index);
			poly_split_fft9(tx, z0, index1, lvl);
			//poly_copy_fft9(ty, z1, lvl, j);
			lvl = lvl - 1;
			idx = idx << 1;
			if (lvl == 1) state = 1;
		}
		if (state == 1) {
			n = 1 << lvl;
			index = convert(idx, 9 - lvl, 1);
			index1 = convert_fft(idx, 9 - lvl);
			i1 = convert3(lvl) + convert2(index, lvl); i2 = i1 + n; i3 = i1 + n + (lvl << (lvl - 1));
			printf("state1 index %d  j %d  index1 %d   lvl %d  i1 %d  i2 %d  i3 %d \n", index, idx, index1, lvl, i1, i2, i3);
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
			printf("state2   index %d  j %d index1 %d  lvl %d  i1 %d  i2 %d  i3 %d conv3 %d  \n", index, idx, index1, lvl, i1, i2, i3, convert3(lvl));
			printf(" \n");
			poly_sub_logn(tx, z0, lvl, index1);
			idx++;
			index = convert(idx, 9 - lvl, 1);
			i1 = convert3(lvl) + convert2(index, lvl); i2 = i1 + n;
			index1 = convert_fft(idx, 9 - lvl);
			for (int i = 0; i < n; i++) {
				treex[i * p + index1] = tree[i1 + i];
				//	printf(" treex[%d]  %lf  ,", i * p + index1, treex[i * p + index1]);
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
			printf("state3  index %d  j %d  index1 %d  lvl %d  i1 %d  i2 %d  i3 %d  \n", index, idx, index1, lvl, i1, i2, i3);
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
