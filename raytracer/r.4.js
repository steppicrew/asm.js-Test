
// FIXME: Uebergeben
var w = 50;

function Raytracer( stdlib, foreign, heap ) {

    // "use asm";

    function run() {

        // Tiny Raytracer (C) Gabriel Gambetta 2013
        // ----------------------------------------
        //
        //  Configuration and scene
        //
        // Size of the canvas. w is also reused as a "big constant" / "+infinity"
        var w = 50;

        w= heap[heap[0]];

        // Sphere: radius, [cx,  cy,  cz], R,  G,  B, specular exponent, reflectiveness 
        // R, G, B in [0, 9], reflectiveness in [0..9].
        var spheres = [
            w,   0, -w, 0,  9, 9, 0,  w,  2,  // Yellow sphere
//            1,   0,  0, 3,  9, 0, 0,  w,  3,  // Red sphere
//            1,  -2,  1, 4,  0, 9, 0,  9,  4,  // Green sphere
//            1,   2,  1, 4,  0, 0, 9,  w,  5   // Blue sphere
        ];

        // Ambient light.
        var ambient_light = 2;

        // Point lights: intensity, [x,  y,  z]
        // Intensities should add to 10, including ambient.
        var lights = [
            8, [2, 2, 0]
        ];

        // -----------------------------------------------------------------------------

        // Shorten some names.
        var math = Math;
        var sqrt = math.sqrt;
        var max = math.max;

        // Closure doesn't rename vars unless they're declared with "var", which takes
        // space. So most vars are 1-letter and global:
        //
        // C: sphere center
        // L: light vector
        // N: surface normal at intersection
        // X: intersection point
        // a: quadratic equation constant
        // b: quadratic equation constant
        // c: color channel
        // d: quadratic equation discriminant
        // e: loop variable
        // f: candidate parameter t
        // h: half-width of the canvas
        // i: illumination
        // J: (ray origin) - (sphere center) 
        // k: <N, L> 
        // l: light index in loop
        // n: <N, N>
        // q: sphere index in loop
        // r: sphere radius
        // s: closest intersection sphere index
        // t: closest intersection t
        // u: intensity of lights[l] 
        // v: closest sphere found in loop
        //
        // The exceptions are vars that need to be initialized here (we still pay the
        // "a=", so we pay a single "var" above, and use nice names) and some vars in 
        // trace_ray, which is recursive, so some of it vars can't be global.

        // Dot product.
        function dot(Ax, Ay, Az, Bx, By, Bz) {
            return Ax * Bx + Ay * By + Az * Bz;
        }


        // Helper: A_minus_Bk(A, B, k)  =  A - B*k. Since it's used more with k < 0,
        // using - here saves a couple of bytes later.
        // function A_minus_Bk (A, B, k) {
        //     return [A[0] - B[0]*k, A[1] - B[1]*k, A[2] - B[2]*k];
        // }

        var t;

        // Find nearest intersection of the ray from B in direction D with any sphere.
        // "Interesting" parameter values must be in the range [t_min, t_max].
        // Returns the index within spheres of the center of the hit sphere, 0 if none.
        // The parameter value for the intersection is in the global variable t.
        function closest_intersection(Bx, By, Bz, Dx, Dy, Dz, t_min, t_max) {
            t = w;  // Min distance found.

            // For each sphere.
            // Get the radius and test for end of array at the same time; 
            // spheres[n] == undefined ends the loop.
            // q points to the 2nd element of the sphere because of q++; +6 skips to next
            // sphere.
            var v, q, r, J, a, b, d, e, f;
            for (v = q = 0; r = spheres[q++]; q += 8) {  
                // Compute quadratic equation coefficients K1, K2, K3
                
                Jx = Bx - spheres[q];      // origin - center
                Jy = By - spheres[q + 1];  // origin - center
                Jz = Bz - spheres[q + 2];  // origin - center
                
                a =  2 * dot(Dx, Dy, Dz, Dx, Dy, Dz);  // 2*K1
                b = -2 * dot(Jx, Jy, Jz, Dx, Dy, Dz);  // -K2
// console.log(r, a, b, Jx, Jy, Jz, Dx, Dy, Dz);

                // Compute sqrt(Discriminant) = sqrt(K2*K2 - 4*K1*K3), go ahead if there are
                // solutions.
                if ( d = sqrt(b * b - 2 * a *(dot(Jx, Jy, Jz, Jx, Jy, Jz) - r * r)) ) {
                    // Compute the two solutions.
                    for (e = 2; e--; d = -d) {
                        f = (b - d)/a;  // f = (-K2 - d) / 2*K1
// console.log("f", b, b * b, 2 * a * dot(Jx, Jy, Jz, Jx, Jy, Jz), r * r, "e", e);
                        if (t_min < f && f < t_max && f < t) { 
                            v = q;
                            t = f;
                        }
                    }
                }
            }

            // Return index of closest sphere in range; t is global
            return v;
        }


        // Trace the ray from B with direction D considering hits in [t_min, t_max].
        // If depth > 0, trace recursive reflection rays.
        // Returns the value of the current color channel as "seen" through the ray.
        function trace_ray(Bx, By, Bz, Dx, Dy, Dz, t_min, t_max, depth) {
// console.log(depth, Dx, Dy, Dz);

            var s, Nx, Ny, Nz, Xx, Xy, Xz, n, i, l, u, k, Lx, Ly, Lz, Mx, My, Mz;

            // Find nearest hit; if no hit, return black.
            if (!(s = closest_intersection(Bx, By, Bz, Dx, Dy, Dz, t_min, t_max))) {
                return 0;
            }
console.log("s", s);

            // Compute "normal" at intersection: N = X - spheres[s]
            // intersection: X = B + D*t = B - D*(-t)
            // X = A_minus_Bk(Bx, By, Bz, Dx, Dy, Dz, -t);
            Xx = Bx + Dx * t;
            Xy = By + Dy * t;
            Xz = Bz + Dz * t;
            
            Nx = Xx - spheres[s];
            Ny = Xy - spheres[s + 1];
            Nz = Xz - spheres[s + 2];

            // Instead of normalizing N, we divide by its length when appropriate. Most of
            // the time N appears twice, so we precompute its squared length.
            n = dot(Nx, Ny, Nz, Nx, Ny, Nz);

            // Start with ambient light only
            i = ambient_light;

            // For each light
            for (l = 0; u = lights[l++]; ) { // Get intensity and check for end of array

                var LIGHT= lights[l++];
            
                // Compute vector from intersection to light (L = lights[l++] - X) and
                // k = <N,L> (reused below)
                // L = A_minus_Bk(lights[l++], Xx, Xy, Xz, 1);
                Lx = LIGHT[0] - Xx;
                Ly = LIGHT[1] - Xy;
                Lz = LIGHT[2] - Xz;

                k = dot(Nx, Ny, Nz, Lx, Ly, Lz);

                // Add to lighting

                // If the point isn't in shadow
                // [t_min, t_max]  =  [epsilon,  1] - epsilon avoids self-shadow, 1 
                // doesn't look farther than the light itself.

                // Diffuse lighting, only if it's facing the point 
                // <N,L> / (|N|*|L|) = cos(alpha)
                // Also, |N|*|L| = sqrt(<N,N>)*sqrt(<L,L>) = sqrt(<N,N>*<L,L>)

                // Specular highlights
                //
                // specular = (<R,V>   / (|R|*|V|))   ^ exponent
                //          = (<-R,-V> / (|-R|*|-V|)) ^ exponent
                //          = (<-R,D>  / (|-R|*|D|))  ^ exponent
                //
                // R = 2*N*<N,L> - L
                // M = -R = -2*N*<N,L> + L = L + N*(-2*<N,L>)
                //
                // If the resultant intensity is negative, treat it as 0 (ignore it).

                var cc= closest_intersection(Xx, Xy, Xz, Lx, Ly, Lz, 1 / w, 1);

                // M = A_minus_Bk(Lx, Ly, Lz, Nx, Ny, Nz, 2*k/n), D) 
                var mf= 2 * k / n;
                Mx = Lx - Nx * mf;
                My = Ly - Ny * mf;
                Mz = Lz - Nz * mf;
                
                i += u * 
                    !cc * (
                        max(0, k / sqrt(dot(Lx, Ly, Lz, Lx, Ly, Lz) * n))
                        + max(0, math.pow(dot(Mx, My, Mz, Dx, Dy, Dz) 
                            / sqrt(dot(Mx, My, Mz, Mx, My, Mz) * dot(Dx, Dy, Dz, Dx, Dy, Dz)), spheres[s + 6]))
                    );
            }


            // Compute the color channel multiplied by the light intensity. 2.8 maps
            // the color range from [0, 9] to [0, 255] and the intensity from [0, 10]
            // to [0, 1],  because 2.8 ~ (255/9)/10
            // 
            // spheres[s] = sphere center, so spheres[s+c] = color channel
            // (c = [1..3] because ++c below)
            var local_color = spheres[s + c + 2] * i * 2.8;

            // If the recursion limit hasn't been hit yet, trace reflection rays.
            // N = normal (non-normalized - two divs by |N| = div by <N,N>
            // D = -view
            // R = 2*N*<N,V>/<N,N> - V = 2*N*<N,-D>/<N,N> + D = D - N*(2*<N,D>/<N,N>)
            var ref = spheres[s + 7] / 9;

            if ( depth-- ) {

                var Rx, Ry, Rz;
                var rf = 2 * dot(Nx, Ny, Nz, Dx, Dy, Dz) / n;
                Rx= Dx - Nx * rf;
                Ry= Dy - Ny * rf;
                Rz= Dz - Nz * rf;

                return trace_ray(Xx, Xy, Xz, Rx, Ry, Rz, 1/w, w, depth) * ref
                    + local_color * (1 - ref);
            }
                    
            return local_color;
        }

        var y, h, x, c;

        var inx_output= heap[4]|0;

        var out_idx = inx_output + 1;

        // For each y; also compute h=w/2 without paying an extra ";"
        for (y = h=w/2; y-- > -h;) {

            // For each x
            for (x = -h; x++ < h;) {

                // One pass per color channel (!). This way we don't have to deal with
                // "colors".
                for (c = 0; ++c < 4;) {
                    // Camera is at (0, 1, 0)
                    //
                    // Ray direction is (x*vw/cw, y*vh/ch, 1) where vw = viewport width, 
                    // cw = canvas width (vh and ch are the same for height). vw is fixed
                    // at 1 so (x/w, y/w, 1)
                    //
                    // [t_min, t_max] = [1, w], 1 starts at the projection plane, w is +inf
                    //
                    // 2 is a good recursion depth to appreciate the reflections without
                    // slowing things down too much
                    //
// window.ct= window.ct ? window.ct + 1 : 1; console.log(window.ct, c, x, y, w, h);
                    heap[out_idx++] = trace_ray(0, 1, 0, x / w, y / w, 1, 1, w, 2);
                }
                heap[out_idx++] = 255; // Opaque alpha
            }
        }

        heap[inx_output]= out_idx - inx_output - 1;
    };

    return {
        run: run
    };
}

function run2_old() {
    var date= new Date();
    console.log("Start");

    var heap= new ArrayBuffer(100000);

    var raytracer= Raytracer(window, {}, heap);
    raytracer.run();

    var raw_data= heap;
    
    // Get to the raw pixel data.
    var canvas = document.getElementById("c");
    var context2d = canvas.getContext("2d");

    canvas.width = canvas.height = w;

    var image_data= context2d.createImageData(w, w);
    var image_data_data= image_data.data;
    var i= heap[0];
    while ( --i >= 0 ) image_data_data[i] = heap[i + 1];

    context2d.putImageData(image_data, 0, 0);

    console.log("End:", new Date() - date);
}

