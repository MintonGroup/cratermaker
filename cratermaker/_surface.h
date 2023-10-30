#define STRMAX 512

typedef struct surface_type {
    double *elevation;
    int elevation_shape[2];
    char stringvar[STRMAX];
}c_surface_type;
extern struct surface_type* bind_surface_init(int ny, int nx);
extern void bind_surface_final(struct surface_type *obj);
extern void bind_surface_set_stringvar(struct surface_type *obj, const char *c_string);
extern char* bind_surface_get_stringvar(struct surface_type *obj);
