//angle=25;
$fn=100;
union() 
{
	difference()
	{
		//translate([250,0,0]) rotate([0,0,0]) cube(700,center=true);
		translate([0, 389/2, 0]) rotate([90, 0, 0]) linear_extrude(height=389) polygon(points=[[-1044,50],[0,50],[0,338-222*sin(angle)],[-222*cos(angle),338],[-1044,338]]);

		//z-plus
		difference() 
		{
			translate([-1394, 0, 688]) rotate([ 0, 0, 0]) cube(900,center=true);
			translate([ -944, 0, 238]) rotate([90, 0, 0]) cylinder(h = 1000, r = 100, center = true);
		}

		//z-minus
		difference() 
		{
			translate([-1394, 0,-300]) rotate([ 0, 0, 0]) cube(900,center=true);
			translate([ -944, 0, 150]) rotate([90, 0, 0]) cylinder(h = 1000, r = 100, center = true);
		}

		//y-minus
		difference() 
		{
			translate([-1394, -544.5, 0]) rotate([0,0,0])     cube(900,center=true);
			translate([ -944,  -94.5, 0]) rotate([0,0,0]) cylinder(h = 1000, r = 100, center = true);
		}

		//y-plus
		difference() 
		{
			translate([-1394, 544.5, 0]) rotate([0, 0, 0])     cube(900, center=true);
			translate([ -944,  94.5, 0]) rotate([0, 0, 0]) cylinder(h = 1000, r = 100, center = true);
		}
	}

	translate([-842,  163.5, 30])  cylinder(h=60, r=15, center=true);
	translate([-842, -163.5, 30])  cylinder(h=60, r=15, center=true);
	translate([-435,  163.5, 30])  cylinder(h=60, r=15, center=true);
	translate([-435, -163.5, 30])  cylinder(h=60, r=15, center=true);
}
