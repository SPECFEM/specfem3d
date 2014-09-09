

#define JB_NL  11
#define JB_NH  15
#define JB_ND  118

typedef struct _JBTable JBTable;
struct _JBTable {
	int    nl;
	float  z[JB_NL];
	float  vp[JB_NL];
	float  vs[JB_NL];
	float  den[JB_NL];

	int   nh;
	float depths[JB_NH];
	int   nd;
	float deltas[JB_ND];
	float ts[JB_NH][JB_ND];
};

float      gfact3(JBTable *jb, float dist, float depth);
JBTable *  jbtable_load();
