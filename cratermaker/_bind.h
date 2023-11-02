#define STRMAX 512

typedef struct c_surface_type {
    double *elevation;
    int elevation_shape[2];
    char stringvar[STRMAX];
}surface_type;
extern struct c_surface_type* bind_surface_init(int ny, int nx);
extern void bind_surface_final(struct c_surface_type *obj);
extern void bind_surface_set_stringvar(struct c_surface_type *obj, const char *c_string);
extern char* bind_surface_get_stringvar(struct c_surface_type *obj);

typedef struct c_PerlinArguments {
      double damp;
      double damp0;
      double damp_scale;
      double freq;
      double gain;
      double gain0;
      double lacunarity;
      double noise_height;
      double pers;
      double slope;
      double warp;
      double warp0;
}PerlinArguments;