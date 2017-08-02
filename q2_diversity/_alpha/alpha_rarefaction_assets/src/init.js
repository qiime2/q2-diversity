/* global categories */
/* global metrics */

import { select } from 'd3';

import { state } from './data';

import { addMetricPicker, addCategoryPicker } from './toolbar';


export default function init() {
  const metric = metrics[0];
  const category = categories[0];
  // DOM
  const body = select('#main');
  const plotRow = body.insert('div', ':first-child').attr('class', 'viz row');
  const plotDiv = plotRow.append('div').attr('class', 'col-lg-12');
  const controlsRow = plotDiv.append('div').attr('class', 'controls row');
  const svgCol = plotRow.append('div').attr('class', 'col-lg-10');
  const svg = svgCol.append('svg');
  const chart = svg.append('g');
  chart.append('g').attr('class', 'x axis');
  chart.append('g').attr('class', 'y axis');
  chart.append('text').attr('class', 'x label');
  chart.append('text').attr('class', 'y label');
  const legend = plotRow.append('div')
                       .attr('class', 'col-lg-2')
                       .style('height', '450px')
                       .style('float', 'left')
                       .style('overflow-y', 'scroll')
                       .style('overflow-x', 'auto')
                       .append('svg')
                       .attr('viewBox', '0 0 200 1000')
                       .attr('class', 'legend')
                       .append('g');
  body.insert('h1', ':first-child').text('Alpha Rarefaction');
  // D3
  state.initialize(metric, category, controlsRow, svg, legend);
  addMetricPicker(controlsRow, metrics, metric);
  if (categories.length > 0) {
    addCategoryPicker(controlsRow, categories, category);
  }
}
