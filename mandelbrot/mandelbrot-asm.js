
function Mandelbrot( stdlib, foreign, heap ) {
    'use asm';

    var data= new stdlib.Uint8Array(heap);
    var imul= stdlib.Math.imul;

    function _data( ctx_tile_x, ctx_tile_y, map_zoom, tileSize ) {
        ctx_tile_x= ctx_tile_x|0;
        ctx_tile_y= ctx_tile_y|0;
        map_zoom= map_zoom|0;
        tileSize= tileSize|0;

        var tileCount = 1;
        var ReStart = -2.0;
        var ReDiff = 3.0;
        var MinRe = 0.0;
        var MaxRe = 0.0;
        var ImStart = -1.2;
        var ImDiff = 2.4;
        var MinIm = 0.0;
        var MaxIm = 0.0;
        var Re_factor = 0.0;
        var Im_factor = 0.0;
        var MaxIterations = 32;
        var MaxIterations2 = 16;

        var i= 0;
        var y= 0;
        var c_im= 0.0;
        var x= 0;
        var c_re= 0.0;
        var Z_re= 0.0;
        var Z_im= 0.0;
        var isInside= 0;
        var n= 0;
        var Z_re2= 0.0;
        var Z_im2= 0.0;

        tileCount = 1 << map_zoom;
        MinRe = +(ReStart + ReDiff * +(ctx_tile_x >> 0) / +(tileCount >> 0));
        MaxRe = MinRe + ReDiff / +(tileCount >> 0);
        MinIm = ImStart + ImDiff * +(ctx_tile_y >> 0) / +(tileCount >> 0);
        MaxIm = MinIm + ImDiff / +(tileCount >> 0);
        Re_factor = (MaxRe - MinRe) / +((tileSize - 1) >> 0);
        Im_factor = (MaxIm - MinIm) / +((tileSize - 1) >> 0);

        i= 0;
        for ( y = 0; (y|0) < (tileSize|0); y= (y + 1)|0 ) {
            c_im = MinIm + (+(y|0)) * Im_factor;
            for ( x = 0; (x|0) < (tileSize|0); x= (x + 1)|0 ) {
                c_re = MinRe + (+(x|0)) * Re_factor;
                Z_re = c_re;
                Z_im = c_im;
                isInside = 1;
                for ( n = 0; (n|0) < (MaxIterations|0); n= (n + 1)|0 ) {
                    Z_re2 = Z_re * Z_re;
                    Z_im2 = Z_im * Z_im;
                    if ( Z_re2 + Z_im2 > +4 ) {
                        isInside = 0;
                        break;
                    }
                    Z_im = +2 * Z_re * Z_im + c_im;
                    Z_re = Z_re2 - Z_im2 + c_re;
                }
                if ( isInside ) {
                    data[i] = data[(i + 1)|0] = data[(i + 2)|0] = 0;
                }
                else if ( (n|0) < (MaxIterations2|0) ) {
                    data[i] = (imul(255, n)|0) / (MaxIterations2|0);
                    data[(i + 1)|0] = data[(i + 2)|0] = 0;
                }
                else {
                    data[i] = 255;
                    data[(i + 1)|0] = data[(i + 2)|0] = (imul((n - MaxIterations2), 255)|0) / (MaxIterations2|0);
                }
                data[(i + 3)|0] = 255;
                i= (i + 4)|0;
            }
        }
    }

    return {
        data: _data,
    };
};



L.MandelbrotSet = L.TileLayer.Canvas.extend({
    tileSize: 256,
    data: undefined,
    mandel: undefined,

    initialize: function (options) {
        L.Util.setOptions(this, options);
        this.drawTile = function (canvas, tilePoint) {
            var ctx = {
                canvas: canvas,
                tile: tilePoint
            };
            this._draw(ctx);
        };
        var data= new ArrayBuffer(this.tileSize * this.tileSize * 4 * 1024);
        this.mandel= Mandelbrot(window, {}, data);
        this.data= new Uint8Array(data);

        var profileTimer;
        this.on('loading', function() {
            profileTimer= new Date();
        })
        this.on('load', function() {
            document.getElementById('profiler').innerHTML= (new Date() - profileTimer) + ' ms<br />';
        })
    },

    _draw: function (ctx) {
        this.mandel.data(ctx.tile.x, ctx.tile.y, this._map._zoom, this.tileSize);
        var data = this.data;
        var g = ctx.canvas.getContext('2d');
        var d = g.createImageData(this.tileSize, this.tileSize);
        for (var i = 0, n = this.tileSize * this.tileSize * 4; i < n; i++)  {
            d.data[i] = data[i];
        }
        g.putImageData(d, 0, 0);
    },
});

var attr = 'Adapted from a <a target="_blank" href="http://polymaps.org/ex/mandelbrot.html">polymaps sample</a>'
var map = L.map('map');
map.attributionControl.addAttribution(attr);
map.setView(new L.LatLng(0, 0), 2);
map.addLayer(new L.MandelbrotSet());

/*
var canvas= document.getElementById('c');
var context2d= canvas.getContext('2d');
canvas.width= 200;
canvas.height= 200;
*/

