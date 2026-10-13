#include <TApplication.h>

void raytracer();

int main(int argc, char **argv) {
  TApplication application("raytracer", &argc, argv);
  raytracer();
  application.Run();
  return 0;
}
