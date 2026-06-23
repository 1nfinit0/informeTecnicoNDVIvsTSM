var distritos = ee.FeatureCollection("users/hsluis4326/distritosPeru"),
  image = ee.Image("users/hsluis4326/logo"),
  icenTable = ee.FeatureCollection("projects/test-hsluis4326/assets/csvTest"),
  modis = ee.ImageCollection("MODIS/061/MOD13Q1"),
  logo = ee.Image("projects/test-hsluis4326/assets/logo"),
  logo2 = ee.Image("projects/test-hsluis4326/assets/logo2");
Map.setOptions("HYBRID");

////////////////////////////////////////////////////////////////////////
// SECCIÓN: PREPROCESAMIENTO DE LA TABLA ICEN Y PROMEDIOS ESTACIONALES
////////////////////////////////////////////////////////////////////////
var withDate = icenTable.map(function (f) {
  var yy = f.get("yy");
  var mm = f.get("mm");
  var numMm = ee.Number.parse(mm);
  var mm2 = ee.Algorithms.If(
    numMm.lt(10),
    ee.String("0").cat(ee.String(numMm)),
    ee.String(numMm),
  );
  var dateStr = ee.String(yy).cat(mm2);
  return f.set("date", dateStr);
});

var sorted = withDate.sort("date", false);
var limitedFc = sorted.limit(12);

var limitedList = limitedFc.toList(12);
var tabla12 = ee.List(
  limitedList.map(function (feat) {
    feat = ee.Feature(feat);
    var yy = ee.String(feat.get("yy"));
    var mm = ee.Number.parse(feat.get("mm"));

    return ee.Dictionary({
      mes: yy.cat("-").cat(mm),
      variabilidad: ee.Number.parse(feat.get("icen")),
    });
  }),
);

var targetYears = ee.List(["2007", "2017", "2023", "2024"]);

var icenPromedios = targetYears.map(function (yy) {
  yy = ee.String(yy);
  var prevY = ee.Number.parse(yy).subtract(1).format();

  var decPrev = icenTable
    .filter(ee.Filter.eq("yy", prevY))
    .filter(ee.Filter.eq("mm", "12"));

  var janToMar = icenTable
    .filter(ee.Filter.eq("yy", yy))
    .filter(ee.Filter.inList("mm", ["1", "2"]));

  var periodFC = decPrev.merge(janToMar);

  var valores = ee.List(periodFC.aggregate_array("icen")).map(function (v) {
    return ee.Number.parse(v);
  });

  var meanICEN = ee.Number(valores.reduce(ee.Reducer.mean()));

  return ee.Dictionary({
    yy: yy,
    icen_promedio: meanICEN,
  });
});

print("Promedio ICEN (dic–mar) para 2007, 2017 y 2024:", icenPromedios);

////////////////////////////////////////////////////////////////////////
// SECCIÓN: PROMEDIO DE ANOMALÍA DE TSM NOAA CDR OISST EN LA REGIÓN NIÑO 1+2
////////////////////////////////////////////////////////////////////////
var regionNino12 = ee.Geometry.Rectangle([-90, -10, -80, 0], null, false);
var anomCol = ee.ImageCollection("NOAA/CDR/OISST/V2_1").select("anom");

var sstAnomPromedios = targetYears.map(function (yyStr) {
  var yearInt = ee.Number.parse(yyStr).toInt();
  var prevYearInt = yearInt.subtract(1).toInt();

  var start = ee.Date.fromYMD(prevYearInt, 12, 1);
  var end = ee.Date.fromYMD(yearInt, 2, 1);

  var anomMeanImg = anomCol.filterDate(start, end).mean().multiply(0.01);

  var dict = anomMeanImg.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: regionNino12,
    scale: 25000,
    maxPixels: 1e13,
  });

  return ee.Dictionary({
    yy: yearInt,
    anom_promedio: dict.get("anom"),
  });
});

print("Promedio anomalía TSM (dic–mar) en Niño 1+2 (°C):", sstAnomPromedios);

////////////////////////////////////////////////////////////////////////
// SECCIÓN: CÁLCULO Y VISUALIZACIÓN DEL NDVI
////////////////////////////////////////////////////////////////////////
var calcularYMostrarNDVIPromedio = function (geometry, yearsBack) {
  var layers = Map.layers();
  for (var i = 0; i < layers.length(); i++) {
    var lyr = layers.get(i);
    if (lyr.getName() === layerName) {
      Map.remove(lyr);
      break;
    }
  }

  var today = ee.Date(Date.now());
  var startDate = today.advance(-yearsBack, "year");

  var s2 = ee
    .ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
    .filterBounds(geometry)
    .filterDate(startDate, today)
    .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", 10))
    .map(function (image) {
      var ndvi = image.normalizedDifference(["B8", "B4"]).rename("NDVI");
      return ndvi;
    });

  var ndvi_mean = s2.mean().clip(geometry);

  var ndviViz = {
    min: -0.2,
    max: 1,
    palette: ["blue", "white", "green"],
  };

  var layerName = "NDVI Promedio (" + yearsBack + " años)";

  Map.addLayer(ndvi_mean, ndviViz, layerName);

  return ndvi_mean;
};

var obtenerAreasAgricolasPermanentes = function (geometry, threshold, scale) {
  var today = ee.Date(Date.now());
  var s2 = ee
    .ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
    .filterBounds(geometry)
    .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", 10))
    .select(["B8", "B4"]);

  var years = ee.List.sequence(0, 3);

  var masks = ee.ImageCollection.fromImages(
    years.map(function (clientOffset) {
      var offset = ee.Number(clientOffset);
      var currentYear = today.get("year");
      var startYear = currentYear.subtract(offset.add(1));
      var start = ee.Date.fromYMD(startYear, 1, 1);
      var end = start.advance(1, "year");

      var annualNDVI = s2
        .filterDate(start, end)
        .map(function (image) {
          return image
            .normalizedDifference(["B8", "B4"])
            .rename("NDVI")
            .clip(geometry);
        })
        .max();

      return annualNDVI.gt(threshold).clip(geometry).rename("mask");
    }),
  );

  var count = masks.reduce(ee.Reducer.sum()).rename("countYears");
  var persistent = count.gte(3).selfMask().clip(geometry);

  var vectors = persistent.reduceToVectors({
    geometry: geometry,
    geometryType: "polygon",
    scale: scale,
    maxPixels: 1e9,
    tileScale: 4,
  });

  Map.addLayer(
    persistent,
    { palette: ["00000000", "00FF00"] },
    "Áreas Agrícolas Permanentes (NDVI > 0.6, 4 años)",
  );

  return vectors;
};

////////////////////////////////////////////////////////////////////////
// SECCIÓN: INTERFAZ DE SELECCIÓN GEOGRÁFICA Y FILTRO TEMPORAL
////////////////////////////////////////////////////////////////////////
var departamentos = distritos.aggregate_array("NOMBDEP").distinct().sort();

var provinciaSelect = ui.Select({
  placeholder: "Selecciona una provincia",
  disabled: true,
});
var distritoSelect = ui.Select({
  placeholder: "Selecciona un distrito",
  disabled: true,
});
var monthMap = {
  1: "Enero",
  2: "Febrero",
  3: "Marzo",
  4: "Abril",
  5: "Mayo",
  6: "Junio",
  7: "Julio",
  8: "Agosto",
  9: "Septiembre",
  10: "Octubre",
  11: "Noviembre",
  12: "Diciembre",
};
var monthItems = Object.keys(monthMap).map(function (code) {
  return { label: monthMap[code], value: code };
});
var mesSelector = ui.Select({
  placeholder: "Selecciona un mes",
  disabled: true,
});

var departamentoSelect = ui.Select({
  items: departamentos.getInfo(),
  placeholder: "Selecciona un departamento",
  onChange: function (departamento) {
    var provincias = distritos
      .filter(ee.Filter.eq("NOMBDEP", departamento))
      .aggregate_array("NOMBPROV")
      .distinct()
      .sort();

    provinciaSelect.items().reset(provincias.getInfo());
    provinciaSelect.setDisabled(false);

    if (!panel.widgets().contains(provinciaSelect)) {
      panel.add(provinciaSelect);
    }

    distritoSelect.setDisabled(true);
    distritoSelect.items().reset([]);
    panel.widgets().remove(distritoSelect);
    mesSelector.setDisabled(true);
    mesSelector.items().reset([]);
    panel.widgets().remove(mesSelector);

    var widgets = panel.widgets();
    var count = widgets.length();
    for (var i = 0; i < count; i++) {
      var widget = widgets.get(i);
      if (
        widget instanceof ui.Label &&
        (widget.getValue() === "Monitoreo ICEN (últimos 12 registros)" ||
          widget.getValue() ===
            "Monitoreo TSM datos de NOAA (últimos 12 meses actualizado)" ||
          widget.getValue() ===
            "NDVI medio vs TSM (registro histórico de 15 años)" ||
          widget.getValue() === "Correlación TSM vs NDVI")
      ) {
        panel.widgets().remove(widget);
        i--;
        count--;
      } else if (widget instanceof ui.Chart) {
        panel.widgets().remove(widget);
        i--;
        count--;
      }
    }
  },
});

provinciaSelect.onChange(function (provincia) {
  var distritosFiltrados = distritos
    .filter(ee.Filter.eq("NOMBPROV", provincia))
    .aggregate_array("NOMBDIST")
    .distinct()
    .sort();

  distritoSelect.items().reset(distritosFiltrados.getInfo());
  distritoSelect.setDisabled(false);

  if (!panel.widgets().contains(distritoSelect)) {
    panel.add(distritoSelect);
  }
});

distritoSelect.onChange(function (distrito) {
  mesSelector.items().reset(monthItems);
  mesSelector.setDisabled(false);
  if (!panel.widgets().contains(mesSelector)) {
    panel.add(mesSelector);
  }

  var distritoSeleccionado = distritos
    .filter(ee.Filter.eq("NOMBDIST", distritoSelect.getValue()))
    .filter(ee.Filter.eq("NOMBPROV", provinciaSelect.getValue()))
    .filter(ee.Filter.eq("NOMBDEP", departamentoSelect.getValue()));

  Map.layers().reset();

  var bordeColor = "red";
  var rellenoColor = "36414900";

  Map.addLayer(
    distritoSeleccionado.style({
      color: bordeColor,
      fillColor: rellenoColor,
      width: 1.5,
    }),
    {},
    "Distrito Seleccionado",
  );

  Map.centerObject(distritoSeleccionado);
});

mesSelector.onChange(function (mes) {
  var widgets = panel.widgets();
  var count = widgets.length();
  for (var i = 0; i < count; i++) {
    var widget = widgets.get(i);
    if (
      widget instanceof ui.Label &&
      (widget.getValue() === "Monitoreo ICEN (últimos 12 registros)" ||
        widget.getValue() ===
          "Monitoreo TSM datos de NOAA (últimos 12 meses actualizado)" ||
        widget.getValue() ===
          "NDVI medio vs TSM (registro histórico de 15 años)" ||
        widget.getValue() === "Correlación TSM vs NDVI")
    ) {
      panel.widgets().remove(widget);
      i--;
      count--;
    } else if (widget instanceof ui.Chart) {
      panel.widgets().remove(widget);
      i--;
      count--;
    }
  }

  var valores12 = tabla12.map(function (d) {
    return ee.Number(ee.Dictionary(d).get("variabilidad"));
  });
  var etiquetas12 = tabla12.map(function (d) {
    return ee.String(ee.Dictionary(d).get("mes"));
  });
  var vals12 = valores12.getInfo();
  var labs12 = etiquetas12.getInfo();
  vals12.reverse();
  labs12.reverse();
  labs12 = labs12.map(function (label) {
    return label.split("-")[1];
  });
  var chart12 = ui.Chart.array
    .values({
      array: vals12,
      axis: 0,
      xLabels: labs12,
    })
    .setSeriesNames(["Variabilidad ICEN"])
    .setChartType("LineChart")
    .setOptions({
      title: "Variabilidad de la TSM datos del ICEN - IGP",
      vAxis: {
        title: "Variación ICEN (°C)",
        viewWindow: { min: -3, max: 3 },
        gridlines: { count: 7 },
        baseline: 0,
        baselineColor: "black",
      },
      hAxis: {
        title: "Mes (De pasado a presente →)",
        gridlines: { count: 12 },
        slantedText: true,
        slantedTextAngle: 45,
      },
      lineWidth: 2,
      pointSize: 4,
      colors: ["steelblue"],
      legend: { position: "none" },
    });

  var titulo12 = ui.Label("Monitoreo ICEN (últimos 12 registros)", {
    fontWeight: "bold",
    fontSize: "14px",
    margin: "20px 0 10px 0",
    textAlign: "center",
  });

  if (panel.widgets().contains(titulo12)) panel.widgets().remove(titulo12);
  if (panel.widgets().contains(chart12)) panel.widgets().remove(chart12);
  panel.add(titulo12);
  panel.add(chart12);

  var valores = tablaAnom.map(function (d) {
    return ee.Number(ee.Dictionary(d).get("variabilidad"));
  });
  var etiquetas = tablaAnom.map(function (d) {
    return ee.String(ee.Dictionary(d).get("mes"));
  });
  var vals = valores.getInfo();
  var labs = etiquetas.getInfo();
  vals.reverse();
  labs.reverse();
  labs = labs.map(function (label) {
    return label.split("-")[1];
  });
  var chartSST = ui.Chart.array
    .values({
      array: vals,
      axis: 0,
      xLabels: labs,
    })
    .setSeriesNames(["Anomalía TSM"])
    .setChartType("LineChart")
    .setOptions({
      title: "Calculado en la zona 1+2",
      vAxis: {
        title: "Anomalía TSM (°C)",
        viewWindow: { min: -3, max: 3 },
        gridlines: { count: 7 },
        baseline: 0,
        baselineColor: "black",
      },
      hAxis: {
        title: "Mes (De pasado a presente ->)",
        gridlines: { count: 12 },
        slantedText: true,
        slantedTextAngle: 45,
      },
      lineWidth: 2,
      pointSize: 4,
      colors: ["steelblue"],
      legend: { position: "none" },
    });

  var tituloSST = ui.Label(
    "Monitoreo TSM datos de NOAA (últimos 12 meses actualizado)",
    {
      fontWeight: "bold",
      fontSize: "14px",
      margin: "20px 0 10px 0",
      textAlign: "center",
    },
  );

  if (panel.widgets().contains(tituloSST)) panel.widgets().remove(tituloSST);
  if (panel.widgets().contains(chartSST)) panel.widgets().remove(chartSST);
  panel.add(tituloSST);
  panel.add(chartSST);

  //////////////////////////////////////////////////////////////////////////
  // SECCIÓN: PARÁMETROS COMUNES DE ANÁLISIS TEMPORAL
  //////////////////////////////////////////////////////////////////////////
  var selMonth = mes;
  var selMonthNum = ee.Number.parse(selMonth);
  var selMonthDate = ee.Date.fromYMD(2000, selMonthNum, 1);
  var selMonthLabel = selMonthDate.format("MM");

  var today = ee.Date(Date.now());
  var currentYear = today.get("year");
  var yearsList = ee.List.sequence(currentYear, currentYear.subtract(14), -1);

  var geom = distritos
    .filter(ee.Filter.eq("NOMBDIST", distritoSelect.getValue()))
    .filter(ee.Filter.eq("NOMBPROV", provinciaSelect.getValue()))
    .filter(ee.Filter.eq("NOMBDEP", departamentoSelect.getValue()));

  var geomBounds = obtenerAreasAgricolasPermanentes(geom, 0.8, 50);
  var geomMulti = geomBounds.geometry();
  var highCount = geomBounds.size();

  //////////////////////////////////////////////////////////////////////////
  // SECCIÓN: SERIE ICEN DESPLAZADA DOS MESES
  //////////////////////////////////////////////////////////////////////////
  var icenLast15 = yearsList.map(function (y) {
    var date = ee.Date.fromYMD(ee.Number(y), selMonthNum, 1);
    var offset = date.advance(-2, "month");
    var yearStr = offset.format("YYYY");
    var monStr = offset.format("M");
    var label = offset.format("YYYY-MM");
    var filtered = icenTable
      .filter(ee.Filter.eq("yy", yearStr))
      .filter(ee.Filter.eq("mm", monStr));
    var raw = ee.Algorithms.If(
      filtered.size().gt(0),
      filtered.first().get("icen"),
      "0",
    );
    return ee.Dictionary({
      mes: label,
      variabilidad: ee.Number.parse(raw),
    });
  });
  //////////////////////////////////////////////////////////////////////////
  // SECCIÓN: SERIE NDVI MODIS PARA EL MES SELECCIONADO
  //////////////////////////////////////////////////////////////////////////
  function generarNDVIPorMes_MODIS(geomBounds, selMonth) {
    var mesNum = ee.Number.parse(selMonth);
    var hoy = ee.Date(Date.now());
    var anioActual = hoy.get("year");
    var mesActual = hoy.get("month");
    var mesExiste = mesNum.lte(mesActual);
    var inicioAnios = ee.Algorithms.If(
      mesExiste,
      anioActual,
      ee.Number(anioActual).subtract(1),
    );
    var listaAnios = ee.List.sequence(
      ee.Number(inicioAnios),
      ee.Number(inicioAnios).subtract(14),
      -1,
    );
    var mesStr = mesNum.format("%02d");

    var lista = listaAnios.map(function (anio) {
      anio = ee.Number(anio);
      var inicio = ee.Date.fromYMD(anio, mesNum, 1);
      var fin = inicio.advance(1, "month");

      var img = ee
        .ImageCollection("MODIS/061/MOD13A3")
        .filterDate(inicio, fin)
        .select("NDVI")
        .first();
      img = ee.Image(ee.Algorithms.If(img, img, ee.Image(0).rename("NDVI")));

      var dict = img.multiply(0.0001).reduceRegion({
        reducer: ee.Reducer.mean(),
        geometry: geomBounds,
        scale: 1000,
        maxPixels: 1e13,
      });
      var ndviMedio = dict.get("NDVI");

      return ee.Dictionary({
        mes: anio.format().cat("-").cat(mesStr),
        ndvi_promedio: ndviMedio,
      });
    });

    var listaConDefault = ee.List(
      ee.Algorithms.If(
        mesExiste,
        lista,
        ee
          .List([
            ee.Dictionary({
              mes: ee.Number(anioActual).format().cat("-").cat(mesStr),
              ndvi_promedio: 0.5,
            }),
          ])
          .cat(lista),
      ),
    );

    var listaCorregida = listaConDefault.map(function (obj) {
      obj = ee.Dictionary(obj);
      var val = obj.get("ndvi_promedio");
      var nuevo = ee.Algorithms.If(ee.Algorithms.IsEqual(val, null), 0.5, val);
      return obj.set("ndvi_promedio", ee.Number(nuevo));
    });

    return ee.List(listaCorregida);
  }

  var promedioNDVIs = generarNDVIPorMes_MODIS(geomBounds, selMonth);
  //////////////////////////////////////////////////////////////////////////
  // SECCIÓN: PROMEDIO NDVI ESTACIONAL POR AÑO
  //////////////////////////////////////////////////////////////////////////
  function promedioNDVIMayJunAug_MODIS(geomBounds) {
    return targetYears.map(function (yyStr) {
      var yearInt = ee.Number.parse(ee.String(yyStr)).toInt();

      var start = ee.Date.fromYMD(yearInt, 5, 1);
      var end = ee.Date.fromYMD(yearInt, 9, 1);

      var imgCol = ee
        .ImageCollection("MODIS/061/MOD13A3")
        .filterDate(start, end)
        .map(function (img) {
          return img.set("month", img.date().get("month"));
        })
        .filter(ee.Filter.inList("month", [5, 6, 8]))
        .select("NDVI")
        .mean()
        .multiply(0.0001);

      var dict = imgCol.reduceRegion({
        reducer: ee.Reducer.mean(),
        geometry: geomBounds,
        scale: 1000,
        maxPixels: 1e13,
      });

      var ndviMedio = ee.Number(dict.get("NDVI", 0));

      return ee.Dictionary({
        yy: yearInt,
        ndvi_promedio: ndviMedio,
      });
    });
  }

  var listaPromedios = promedioNDVIMayJunAug_MODIS(geomBounds);
  print(
    "NDVI promedio (May, Jun, Aug) para 2007, 2017, 2023, 2024:",
    listaPromedios,
  );

  //////////////////////////////////////////////////////////////////////////
  // SECCIÓN: ANOMALÍA SST NOAA DESPLAZADA DOS MESES
  //////////////////////////////////////////////////////////////////////////
  var anomLast15 = yearsList.map(function (y) {
    var yearNum = ee.Number(y);
    var start = ee.Date.fromYMD(yearNum.subtract(1), 12, 1);
    var end = ee.Date.fromYMD(yearNum, 2, 1);

    var stat = ee
      .ImageCollection("NOAA/CDR/OISST/V2_1")
      .select("anom")
      .filterBounds(regionNino12)
      .filterDate(start, end)
      .mean()
      .reduceRegion({
        reducer: ee.Reducer.mean(),
        geometry: regionNino12,
        scale: 25000,
        maxPixels: 1e8,
      });

    var variab = ee.Number(stat.get("anom", 0)).multiply(0.01);
    var label = ee
      .String(start.format("YYYY-MM"))
      .cat("_")
      .cat(ee.Date.fromYMD(yearNum, 1, 1).format("YYYY-MM"));

    return ee.Dictionary({
      mes: label,
      variabilidad: variab,
    });
  });

  print("NOAA TSM (Dic prev + Ene actual):", anomLast15);

  //////////////////////////////////////////////////////////////////////////
  // SECCIÓN: COMPARACIÓN NDVI VS SST Y CORRELACIÓN
  //////////////////////////////////////////////////////////////////////////
  var zipped = ee.List(promedioNDVIs).zip(ee.List(anomLast15));
  var pairedFC = ee.FeatureCollection(
    zipped.map(function (pair) {
      pair = ee.List(pair);
      var nd = ee.Dictionary(pair.get(0));
      var an = ee.Dictionary(pair.get(1));
      return ee.Feature(null, {
        mes: nd.get("mes"),
        ndvi_promedio: nd.get("ndvi_promedio"),
        variabilidad: an.get("variabilidad"),
      });
    }),
  );

  var comboChart = ui.Chart.feature
    .byFeature(pairedFC.sort("mes"), "mes", ["ndvi_promedio", "variabilidad"])
    .setChartType("ComboChart")
    .setOptions({
      title: "NDVI vs anomalía TSM (dos meses de diferencia)",
      hAxis: {
        title: "Muestra mensual",
        slantedText: true,
        slantedTextAngle: 45,
        textStyle: { fontSize: 8 },
      },
      series: {
        0: { type: "bars", targetAxisIndex: 0 },
        1: { type: "line", targetAxisIndex: 1, pointSize: 4 },
      },
      vAxes: {
        0: { title: "NDVI promedio", viewWindow: { min: 0.3, max: 1 } },
        1: {
          title: "Anomalía TSM (°C)",
          viewWindow: { min: -3, max: 3 },
          baseline: 0,
        },
      },
      legend: { position: "right" },
      bar: { groupWidth: "90%" },
      chartArea: { width: "80%", height: "70%" },
    });

  var titulo = ui.Label("NDVI medio vs TSM (registro histórico de 15 años)", {
    fontWeight: "bold",
    fontSize: "14px",
  });
  if (panel.widgets().contains(titulo)) panel.widgets().remove(titulo);
  if (panel.widgets().contains(comboChart)) panel.widgets().remove(comboChart);
  panel.add(titulo);
  panel.add(comboChart);

  var pairedList = ee
    .List(promedioNDVIs)
    .zip(ee.List(anomLast15))
    .map(function (pair) {
      pair = ee.List(pair);
      var nd = ee.Dictionary(pair.get(0));
      var an = ee.Dictionary(pair.get(1));
      return ee.Dictionary({
        mes: nd.get("mes"),
        ndvi_promedio: nd.get("ndvi_promedio"),
        variabilidad: an.get("variabilidad"),
      });
    });

  var pairedIcenNdvi = ee
    .List(sstAnomPromedios)
    .zip(ee.List(listaPromedios))
    .map(function (pair) {
      pair = ee.List(pair);
      var ic = ee.Dictionary(pair.get(0));
      var nd = ee.Dictionary(pair.get(1));
      return ee.Feature(null, {
        yy: ic.get("yy"),
        anom_promedio: ic.get("anom_promedio"),
        ndvi_promedio: nd.get("ndvi_promedio"),
      });
    });
  var fc = ee.FeatureCollection(pairedIcenNdvi);

  var corr = fc
    .reduceColumns({
      reducer: ee.Reducer.pearsonsCorrelation(),
      selectors: ["anom_promedio", "ndvi_promedio"],
    })
    .get("correlation");

  var scatterWithTrend = ui.Chart.feature
    .byFeature(fc, "anom_promedio", ["ndvi_promedio"])
    .setChartType("ComboChart")
    .setOptions({
      title:
        "Coeficiente de correlación (r): " +
        ee.Number(corr).format("%.2f").getInfo(),
      hAxis: { title: "TSM promedio (°C)" },
      vAxis: { title: "NDVI promedio" },
      series: {
        0: { type: "scatter" },
      },
      trendlines: {
        0: {
          type: "linear",
          visibleInLegend: false,
          showR2: true,
          lineWidth: 2,
        },
      },
      annotations: {
        textStyle: { fontSize: 12 },
      },
      chartArea: { width: "70%", height: "70%" },
      legend: { position: "none" },
    });

  if (panel.widgets().contains(chartTitle)) panel.remove(chartTitle);
  if (panel.widgets().contains(scatterWithTrend))
    panel.remove(scatterWithTrend);

  var chartTitle = ui.Label("Correlación TSM vs NDVI", {
    fontWeight: "bold",
    fontSize: "14px",
    margin: "0 0 6px 0",
  });

  panel.add(chartTitle);
  panel.add(scatterWithTrend);

  var pairedNdviSst = ee
    .List(listaPromedios)
    .zip(ee.List(sstAnomPromedios))
    .map(function (pair) {
      pair = ee.List(pair);
      var nd = ee.Dictionary(pair.get(0));
      var st = ee.Dictionary(pair.get(1));
      return ee.Dictionary({
        yy: nd.get("yy"),
        ndvi_promedio: nd.get("ndvi_promedio"),
        anom_promedio: st.get("anom_promedio"),
      });
    });

  var fc2 = ee.FeatureCollection(
    pairedNdviSst.map(function (d) {
      return ee.Feature(null, ee.Dictionary(d));
    }),
  );

  var corr2 = fc2
    .reduceColumns({
      reducer: ee.Reducer.pearsonsCorrelation(),
      selectors: ["ndvi_promedio", "anom_promedio"],
    })
    .get("correlation");

  print("Coeficiente de correlación Pearson (NDVI vs anomalía SST)", corr2);
});

var regionNino12 = ee.Geometry.Rectangle({
  coords: [-90, -10, -80, 0],
  proj: "EPSG:4326",
  geodesic: false,
});

var hoy = ee.Date(Date.now());
var meses = ee.List.sequence(0, 11);

var getMonthlyRecord = function (n) {
  n = ee.Number(n);
  var mesDate = hoy.advance(n.multiply(-1), "month");
  var inicio = ee.Date(mesDate.format("YYYY-MM").cat("-01"));
  var fin = inicio.advance(1, "month");

  var col = ee
    .ImageCollection("NOAA/CDR/OISST/V2_1")
    .select("anom")
    .filterBounds(regionNino12)
    .filterDate(inicio, fin);

  var anomMeanImg = col.mean();

  var red = anomMeanImg.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: regionNino12,
    scale: 25000,
    maxPixels: 1e8,
  });

  var raw = ee.Algorithms.If(red.contains("anom"), red.get("anom"), 0);
  var variab = ee.Number(raw).multiply(0.01);

  var etiqueta = inicio.format("YYYY-MM");

  return ee.Dictionary({
    mes: etiqueta,
    variabilidad: variab,
  });
};

var tablaAnom = meses.map(getMonthlyRecord);

var logoThumb = ui.Thumbnail({
  image: logo2.visualize({
    bands: ["b1", "b2", "b3"],
    min: 0,
    max: 255,
  }),
  params: {
    dimensions: "546x121",
    format: "png",
  },
  style: {
    width: "100%",
    height: "auto",
    padding: "0",
    margin: "10px 0",
  },
});

var titulo = ui.Label(
  "MONITOREO DE LA VARIACIÓN EN LA TEMPERATURA SUPERFICIAL DEL MAR DE LA COSTA PERUANA Y LA COBERTURA AGRÍCOLA PERMANENTE POR DISTRITOS",
  {
    fontWeight: "bold",
    fontSize: "16px",
    textAlign: "center",
    margin: "0 auto",
  },
);

var descripcion = ui.Label({
  value:
    "Este análisis evalúa el daño en áreas agrícolas basado en las variaciones de NDVI a lo largo del tiempo. Utilizando imágenes satelitales, se identifican cambios significativos en la vegetación que se correlacionan con las anomalías en la TSM.",
  style: {
    fontSize: "12px",
    textAlign: "center",
    margin: "0 auto",
    padding: "20px",
    width: "350px",
  },
});

var textoEnlace = ui.Label({
  value: "Puede descargar aquí más información",
  targetUrl:
    "https://raw.githubusercontent.com/1nfinit0/informeTecnicoNDVIvsTSM/refs/heads/master/Ficha_Técnica_2025-GEE.pdf",
  style: {
    fontSize: "12px",
    textAlign: "center",
    margin: "0 auto",
    padding: "20px",
    width: "350px",
    color: "blue",
    textDecoration: "underline",
  },
});

var zonaLabel = ui.Label("Zona de interés:", {
  fontWeight: "bold",
  fontSize: "14px",
});

var panel = ui.Panel({
  widgets: [
    logoThumb,
    titulo,
    descripcion,
    textoEnlace,
    zonaLabel,
    departamentoSelect,
  ],
  layout: ui.Panel.Layout.flow("vertical"),
  style: {
    width: "450px",
    padding: "10px",
    margin: "10px auto",
  },
});

ui.root.add(panel);
