/* global d */
import state from './state';

export function addMetricPicker(row, metrics, selectedMetric) {
  row.select('.metricPicker select')
    .on('change', function changeColumn() {
      const newMetric = metrics[this.selectedIndex];
      state.setMetric(newMetric);
    })
    .selectAll('option')
    .data(metrics)
    .enter()
      .append('option')
      .attr('value', d => d)
      .text(d => d)
      .property('selected', d => (d === selectedMetric));
}

export function addColumnPicker(row, columns, selectedColumn) {
  row.select('.columnPicker select')
    .on('change', function changeColumn() {
      const newColumn = columns[this.selectedIndex];
      state.setColumn(newColumn);
    })
    .selectAll('option')
    .data(columns)
    .enter()
      .append('option')
      .attr('value', d => d)
      .text(d => d)
      .property('selected', d => (d === selectedColumn));
}
