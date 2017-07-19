/* global categories */
/* global metrics */

import { select } from 'd3';

import { state } from './data';

import { addMetricPicker, addCategoryPicker } from './toolbar';


export default function init() {
  const metric = metrics[0];
  let category = categories[0];
  if (!category) {
    category = '';
  }
  // DOM
  const body = select('#main');
  const plotRow = body.insert('div', ':first-child').attr('class', 'viz row');
  const plotDiv = plotRow.append('div').attr('class', 'col-lg-12');
  const controlsRow = plotDiv.append('div').attr('class', 'controls row');
  const svgRow = plotDiv.append('div').attr('class', 'plot row');
  const svgCol = svgRow.append('div').attr('class', 'col-lg-12');
  const svg = svgCol.append('svg');
  const chart = svg.append('g');
  body.insert('h1', ':first-child').text('Alpha Rarefaction');
  chart.append('g').attr('class', 'x axis');
  chart.append('g').attr('class', 'y axis');
  chart.append('text').attr('class', 'x label');
  chart.append('text').attr('class', 'y label');

  // D3
  state.initialize(metric, category, controlsRow, svg);
  addMetricPicker(controlsRow, metrics, metric);
  if (categories.length > 0) {
    addCategoryPicker(controlsRow, categories, category);
  }

  // STATS
  // stats(body, data);

  // WARNINGS
  // warnings(body, data);
}
