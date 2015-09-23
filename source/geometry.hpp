#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP 1

class Geometry
{
public:

  virtual double calcVolume(double radius) const = 0;

  virtual double calcArea(double radius) const = 0;

  virtual ~Geometry(void);
};

class Planar: public Geometry
{
public:

  Planar(void);

  double calcVolume(double radius) const;

  double calcArea(double radius) const;
};

class Cylindrical: public Geometry
{
public:

  Cylindrical(void);

  double calcVolume(double radius) const;

  double calcArea(double radius) const;
};

class Spherical: public Geometry
{
public:

  Spherical(void);

  double calcVolume(double radius) const;

  double calcArea(double radius) const;
};

#endif // GEOMETRY_HPP
